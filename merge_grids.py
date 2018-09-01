#!/usr/bin/env python

import argparse
import itertools
import os
import re
import netCDF4
import numpy
import logging

log = logging.getLogger(__name__)


def find_files(data_dir, regex_str, axes):
    regex = re.compile(regex_str)
    file_mapping = {}
    for filepath in os.listdir(data_dir):
        if re.match(regex, os.path.basename(filepath)):
            result = re.search(regex, os.path.basename(filepath))
            key = tuple([result.group(a) for a in axes])
            file_mapping[key] = os.path.join(data_dir, filepath)
    found_keys = file_mapping.keys()
    key_values = [sorted(list(set([k[i] for k in found_keys]))) for i in range(len(axes))]
    for key in itertools.product(*key_values):
        if key not in found_keys:
            numbers = tuple([str(key[axes[i] - 1]).zfill(3) for i in range(len(axes))])
            log.error("Keys not found: %s" % str(numbers))
            return [], {}
    print "The key values are: ", key_values
    return key_values, file_mapping


# TODO: if no experiment number given, just look for all...
def find_xsec(data_dir, expnr, s1, s2, basename):
    regex2d = re.compile("^" + basename + s1 + s2 + ".x(\d+)y(\d+)." + str(expnr).zfill(3) + ".nc$")
    regex3d = re.compile("^" + basename + s1 + s2 + ".(\d+).x(\d+)y(\d+)." + str(expnr).zfill(3) + ".nc$")
    file_mapping_2d, file_mapping_3d = {}, {}
    for filepath in os.listdir(data_dir):
        if re.match(regex2d, os.path.basename(filepath)):
            result = re.search(regex2d, os.path.basename(filepath))
            key = (int(result.group(1)), int(result.group(2)))
            file_mapping_2d[key] = os.path.join(data_dir, filepath)
        elif re.match(regex3d, os.path.basename(filepath)):
            result = re.search(regex3d, os.path.basename(filepath))
            key = (int(result.group(2)), int(result.group(3)), int(result.group(1)))
            file_mapping_3d[key] = os.path.join(data_dir, filepath)
    keys = file_mapping_3d.keys()
    if any(keys):
        dims1, dims2 = max([k[0] for k in keys]) + 1, max([k[1] for k in keys]) + 1
        levs = sorted(list(set(k[2] for k in keys)))
        for i in range(dims1):
            for j in range(dims2):
                for l in levs:
                    if (i, j, l) not in keys:
                        f = '.'.join([basename + s1 + s2, str(l).zfill(4), str(i).zfill(3) + str(j).zfill(3),
                                      str(expnr).zfill(3), "nc"])
                        raise Exception("Missing file: %s" % f)
        return range(dims1), range(dims2), levs, file_mapping_3d
    keys = file_mapping_2d.keys()
    if any(keys):
        dims1, dims2 = max([k[0] for k in keys]) + 1, max([k[1] for k in keys]) + 1
        for i in range(dims1):
            for j in range(dims2):
                if (i, j) not in keys:
                    f = '.'.join([basename + s1 + s2, str(i).zfill(3) + str(j).zfill(3), str(expnr).zfill(3), "nc"])
                    raise Exception("Missing file: %s" % f)
        return range(dims1), range(dims2), [], file_mapping_2d
    return [], [], [], {}


def match_dim_character(varname, ncvar, s):
    index = -1
    counter = 0
    for dim in ncvar.dimensions:
        if dim.startswith(s):
            if index != -1:
                raise Exception("Multiple %s-dimensions found for variable %s" % (s, varname))
            index = counter
        counter += 1
    return index


def build_xsec(dims1, dims2, level_list, file_mapping):
    dst = netCDF4.Dataset("out.nc", 'w')

    # Copy attributes
    datasets, src = {}, None
    for k in file_mapping.keys():
        datasets[k] = netCDF4.Dataset(file_mapping[k], 'r')
        if src is None:
            src = datasets[k]
    dst.setncatts(src.__dict__)

    # Copy and extend spatial dimensions
    dst_dims = {s: 0 for s in src.dimensions.keys() if not src.dimensions[s].isunlimited()}
    for name in dst_dims.keys():
        if not name.startswith('x') and not name.startswith('y'):
            dst_dims[name] = len(src.dimensions[name])
    for i in dims1:
        for j in dims2:
            key = (i, j) if not any(level_list) else (i, j, level_list[0])
            ds = datasets[key]
            for name, dim in ds.dimensions.items():
                if (i == 0 and name.startswith('x')) or (j == 0 and name.startswith('y')):
                    dst_dims[name] += len(dim)
    for name, dim in src.dimensions.items():
        dst.createDimension(name, (dst_dims[name] if not dim.isunlimited() else None))
    if any(level_list):
        levdim = dst.createDimension("lev", len(level_list))
        dst_dims["lev"] = len(levdim)

    # num_steps = len(src.dimensions["time"])
    num_steps = 100
    dt = 10

    # Copy variables
    if any(level_list):
        levvar = dst.createVariable("lev", numpy.float64, ("lev"))
        levvar[:] = numpy.array(level_list)
    dst_vars, time_indices, time_slices = {}, {}, {}
    for name, variable in src.variables.items():
        if any(level_list) and variable.dimensions != (name,):  # Skip axes
            if variable.dimensions[0] == "time":
                dims = ("time", "lev") + variable.dimensions[1:]
            else:
                dims = ("lev",) + variable.dimensions[:]
        else:
            dims = variable.dimensions
        dst_vars[name] = dst.createVariable(name, variable.datatype, dims)
        dst_vars[name].setncatts(src.variables[name].__dict__)
        time_indices[name] = match_dim_character(name, variable, "time")
        if time_indices[name] > 0:
            log.error("Variable %s has non-major time index at %d... skipping variable" % (name, time_indices[name]))
            continue
        if name != "time":
            shape = [dst_dims[d] for d in dims if d in dst_dims]
            if time_indices[name] == 0 and dt > 1:
                if dims[1] == "lev":
                    shape.insert(1, dt)
                else:
                    shape.insert(0, dt)
            time_slices[name] = numpy.full(shape=shape, dtype=numpy.float64, fill_value=numpy.NaN)

    for i in range(0, num_steps, dt):  # loop over t
        istart, iend = i, min(i + dt, num_steps)
        dst_vars["time"][i:iend] = src.variables["time"][i:iend]
        log.info("processing time step %d of %d...", i, num_steps)
        levs = [-1] if not any(level_list) else level_list
        for lev in levs:  # loop over levels
            for j in dims2:  # loop over y
                for k in dims1:  # loop over x
                    key = (k, j) if not any(level_list) else (k, j, lev)
                    ds = datasets[key]
                    for varname, vardata in ds.variables.items():
                        if varname not in time_slices or (time_indices[varname] < 0 and i > 0):
                            continue
                        if time_indices[varname] == 0 and len(ds.dimensions["time"]) == 0:
                            continue
                        time_index = time_indices[varname]
                        axes = (match_dim_character(varname, vardata, 'x'), match_dim_character(varname, vardata, 'y'))
                        if time_index >= 0:
                            values = vardata[...].take(indices=tuple(range(istart, iend)), axis=time_index)
                            if dt == 1:
                                axes = (axes[0] - 1 if time_index < axes[0] else axes[0],
                                        axes[1] - 1 if time_index < axes[1] else axes[1])
                        else:
                            values = vardata[...]
                        if any(level_list) and "lev" in dst_vars[varname].dimensions:
                            zaxis = 0
                            axes = (axes[0] + 1 if axes[0] >= 0 else axes[0], axes[1] + 1 if axes[1] else axes[1])
                            multi_index = (key[0], key[1], level_list.index(lev))
                        else:
                            zaxis = -1
                            multi_index = key
                        copy_block(time_slices[varname], values, multi_index, axes, zaxis)
        for varname, vardata in time_slices.items():
            if time_indices[varname] < 0 and i > 0:
                continue
            if time_indices[varname] == 0:
                if any(level_list):
                    dst_vars[varname][istart:iend, :] = numpy.swapaxes(time_slices[varname], 0, 1)
                else:
                    dst_vars[varname][istart:iend, :] = time_slices[varname]
            elif time_indices[varname] < 0:
                dst_vars[varname][:] = time_slices[varname]
    dst.close()


def copy_block(dest, src, key, xy_axes, z_axis):
    slices = []
    for i in range(len(dest.shape)):
        if i == xy_axes[0]:
            j = i - 1 if 0 <= z_axis < i else i
            slices.append(slice(key[0] * src.shape[j], (key[0] + 1) * src.shape[j], 1))
        elif i == xy_axes[1]:
            j = i - 1 if 0 <= z_axis < i else i
            slices.append(slice(key[1] * src.shape[j], (key[1] + 1) * src.shape[j], 1))
        elif i == z_axis:
            slices.append(slice(key[2], key[2] + 1, 1))
        else:
            slices.append(slice(0, src.shape[i], 1))
    dest[tuple(slices)] = src[...]


def main():
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("datadir", metavar="DIR", type=str, help="Dales output (run) directory")
    parser.add_argument("--exp", "-e", metavar="N", type=int, default=1, help="Experiment number (default: 001)")
    args = parser.parse_args()
    data_dir = args.datadir
    #    dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'x', 'y', "cross")
    #    build_xsec(dims1, dims2, levs, files)
    #    dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'y', 'z', "cross")
    #    build_xsec(dims1, dims2, levs, files)
    #    dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'x', 'z', "cross")
    #    build_xsec(dims1, dims2, levs, files)
    # dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'x', 'y', "surf_")
    # build_xsec(dims1, dims2, levs, files)

    print "Searching 2d cross sections"
    regex = "^crossxy.x(\d+)y(\d+)." + str(args.exp).zfill(3) + ".nc$"
    find_files(data_dir, regex, axes=(1, 2))
    print "Searching 3d cross sections"
    regex = "^crossxy.(\d+).x(\d+)y(\d+)." + str(args.exp).zfill(3) + ".nc$"
    find_files(data_dir, regex, axes=(3, 1, 2))
    print "Searching surface fields"
    regex = "^surf_xy.x(\d+)y(\d+)." + str(args.exp).zfill(3) + ".nc$"
    find_files(data_dir, regex, axes=(1, 2))
    print "Searching 3d fields"
    regex = "^fielddump.(\d+).(\d+)." + str(args.exp).zfill(3) + ".nc$"
    find_files(data_dir, regex, axes=(1, 2))


logging.basicConfig(level=logging.DEBUG)

if __name__ == "__main__":
    main()
