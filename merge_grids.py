#!/usr/bin/env python

import argparse
import os
import re
import netCDF4
import numpy
import logging

log = logging.getLogger(__name__)


# TODO: if no experiment number given, just look for all...
def find_xsec(data_dir, expnr, s1, s2):
    regex2d = re.compile("^cross" + s1 + s2 + ".x(\d+)y(\d+)." + str(expnr).zfill(3) + ".nc$")
    regex3d = re.compile("^cross" + s1 + s2 + ".(\d+).x(\d+)y(\d+)." + str(expnr).zfill(3) + ".nc$")
    file_mapping_2d, file_mapping_3d = {}, {}
    for filepath in os.listdir(data_dir):
        if re.match(regex2d, os.path.basename(filepath)):
            result = re.search(regex3d, os.path.basename(filepath))
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
                        f = '.'.join(["cross" + s1 + s2, str(l).zfill(4), str(i).zfill(3) + str(j).zfill(3),
                                      str(expnr).zfill(3), "nc"])
                        raise Exception("Missing file: %s" % f)
        return range(dims1), range(dims2), levs, file_mapping_3d
    keys = file_mapping_2d.keys()
    if any(keys):
        dims1, dims2 = max([k[0] for k in keys]) + 1, max([k[1] for k in keys]) + 1
        for i in range(dims1):
            for j in range(dims2):
                if (i, j) not in keys:
                    f = '.'.join(["cross" + s1 + s2, str(i).zfill(3) + str(j).zfill(3), str(expnr).zfill(3), "nc"])
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
    dst = netCDF4.Dataset("new.nc", 'w')

    # Copy attributes
    datasets, src = {}, None
    for k in file_mapping.keys():
        datasets[k] = netCDF4.Dataset(file_mapping[k], 'r')
        if src is None:
            src = datasets[k]
    dst.setncatts(src.__dict__)

    # Copy and extend spatial dimensions
    dst_dims = {s: 0 for s in src.dimensions.keys() if not src.dimensions[s].isunlimited()}
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
        dst.createDimension('z', len(level_list))

    # Copy variables
    dst_vars, time_indices, time_slices = {}, {}, {}
    for name, variable in src.variables.items():
        dst_vars[name] = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst_vars[name].setncatts(src.variables[name].__dict__)
        time_indices[name] = match_dim_character(name, variable, "time")
        if time_indices[name] > 0:
            log.error("Variable %s has non-major time index at %d... skipping variable" % (name, time_indices[name]))
            continue
        if name != "time":
            shape = [dst_dims[d] for d in variable.dimensions if d in dst_dims]
            time_slices[name] = numpy.full(shape=shape, dtype=numpy.float64, fill_value=numpy.NaN)

    #    num_steps = len(src.dimensions["time"])
    num_steps = 20
    chunk_size = 10

    for i in range(num_steps):  # loop over t
        dst_vars["time"][i] = src.variables["time"][i]
        print "processing time step", i, "of", num_steps
        for j in dims2:  # loop over y
            for k in dims1:  # loop over x
                if k == 0:
                    continue
                key = (k, j) if not any(level_list) else (k, j, level_list[0])
                ds = datasets[key]
                for varname, vardata in ds.variables.items():
                    if varname not in time_slices or (time_indices[varname] < 0 and i > 0):
                        continue
                    axes = (match_dim_character(varname, vardata, 'x'), match_dim_character(varname, vardata, 'y'))
                    if time_indices[varname] >= 0:
                        values = numpy.asarray(vardata).take(indices=i, axis=time_indices[varname])
                    else:
                        values = vardata[...]
                    copy_xy_block(time_slices[varname], values, key, axes)
        for varname, vardata in time_slices.items():
            if time_indices[varname] < 0 and i > 0:
                continue
            print varname, dst_vars[varname].shape, time_slices[varname].shape
            if time_indices[varname] == 0:
                dst_vars[varname][i, :] = time_slices[varname]
            elif time_indices[varname] < 0:
                dst_vars[varname][:] = time_slices[varname]
        if i % chunk_size == 0:
            dst.sync()
    dst.close()


def copy_xy_block(dest, src, key, axes):
    slices = []
    for i in range(len(dest.shape)):
        found_in_axes = False
        for j in range(len(axes)):
            if i == axes[j]:
                slices.append(slice(start=key[j] * src.shape[j], stop=(key[j] + 1) * src.shape[j]))
                found_in_axes = True
        if not found_in_axes:
            slices.append(slice(start=0, stop=dest.shape[i]))
    dest[tuple(slices)] = src[...]


def main():
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("datadir", metavar="DIR", type=str, help="Dales output (run) directory")
    parser.add_argument("--exp", "-e", metavar="N", type=int, default=1, help="Experiment number (default: 001)")
    args = parser.parse_args()
    data_dir = args.datadir
    dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'x', 'y')
    build_xsec(dims1, dims2, levs, files)


if __name__ == "__main__":
    main()
