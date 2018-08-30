#!/usr/bin/env python

import argparse
import os
import re
import netCDF4
import numpy


#TODO: if no experiment number given, just look for all...
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
    print dst_dims
    for i in dims1:
        for j in dims2:
            key = (i, j) if not any(level_list) else (i, j, level_list[0])
            ds = datasets[key]
            for name, dim in ds.dimensions.items():
                if (i == 0 and name.startswith('x')) or (j == 0 and name.startswith('y')):
                    dst_dims[name] += len(dim)
#                elif name in dst_dims:
#                    if dst_dims[name] != len(dim):
#                        raise Exception("File block dimension %s does not fit lentgh %d" % (name, len(dim)))
    for name, dim in src.dimensions.items():
        dst.createDimension(name, (dst_dims[name] if not dim.isunlimited() else None))
    if any(level_list):
        dst.createDimension('z', len(level_list))

    # Copy variables
    dst_vars, time_indices = {}, {}
    for name, variable in src.variables.items():
        dst_vars[name] = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst_vars[name].setncatts(src.variables[name].__dict__)
        time_indices[name] = match_dim_character(name, variable, "time")



#    num_steps = len(src.dimensions["time"])
    num_steps = 1000
    chunk_size = 10



    for i in range(num_steps):  # loop over t
        print "processing time step", i, "of", num_steps
        var_blocks = {}
        for j in dims2:  # loop over y
            var_slices = {}
            key = (0, j) if not any(level_list) else (0, j, level_list[0])
            ds = datasets[key]
            for varname, vardata in ds.variables.items():
                if varname == "time":
                    continue
                if time_indices[varname] < 0:
                    if i == 0:
                        var_slices[varname] = vardata
                else:
                    var_slices[varname] = numpy.asarray(vardata).take(indices=i, axis=time_indices[varname])
            for k in dims1:  # loop over x
                if k == 0:
                    continue
                key = (k, j) if not any(level_list) else (k, j, level_list[0])
                ds = datasets[key]
                for varname, vardata in ds.variables.items():
                    if varname == "time" or (time_indices[varname] < 0 and i > 0):
                        continue
                    axis_index = match_dim_character(varname, vardata, 'x')
                    if axis_index >= 0:
                        if time_indices[varname] >= 0:
                            values = numpy.asarray(vardata).take(indices=i, axis=time_indices[varname])
                        else:
                            values = numpy.asarray(vardata)
                        axis = axis_index if axis_index < time_indices[varname] else axis_index - 1
                        var_slices[varname] = numpy.append(var_slices[varname], values, axis=axis)
            for varname, vardata in var_slices.items():
                if varname == "time" or (time_indices[varname] < 0 and i > 0):
                    continue
                if varname not in var_blocks:
                    var_blocks[varname] = vardata
                else:
                    axis_index = match_dim_character(varname, src.variables[varname], 'y')
                    if axis_index >= 0:
                        axis = axis_index if axis_index < time_indices[varname] else axis_index - 1
                        var_blocks[varname] = numpy.append(var_blocks[varname], var_slices[varname], axis=axis)
        for varname, vardata in var_blocks.items():
            if time_indices[varname] < 0 and i > 0:
                continue
            print varname, dst_vars[varname].shape, var_blocks[varname].shape
            if time_indices[varname] == 0:
                dst_vars[varname][i, :] = var_blocks[varname]
            else:
                dst_vars[varname][:] = var_blocks[varname]
        if i % chunk_size == 0:
            dst.sync()
    dst.close()


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
