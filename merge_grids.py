#!/usr/bin/env python

import argparse
import os
import re
import netCDF4
import numpy


def find_xsec(data_dir, expnr, s1, s2):
    regex2d = re.compile("^cross" + s1 + s2 + "." + s1 + "(\d+)" + s2 + "(\d+)." + str(expnr).zfill(3) + ".nc$")
    regex3d = re.compile("^cross" + s1 + s2 + ".(\d+)." + s1 + "(\d+)" + s2 + "(\d+)." + str(expnr).zfill(3) + ".nc$")
    file_mapping_2d, file_mapping_3d = {}, {}
    for filepath in os.listdir(data_dir):
        if re.match(regex2d, os.path.basename(filepath)):
            result = re.search(regex3d, os.path.basename(filepath))
            key = (int(result.group(1)), int(result.group(2)))
            file_mapping_2d[key] = filepath
        elif re.match(regex3d, os.path.basename(filepath)):
            result = re.search(regex3d, os.path.basename(filepath))
            key = (int(result.group(2)), int(result.group(3)), int(result.group(1)))
            file_mapping_3d[key] = filepath
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


def build_xsec(s1, dims1, s2, dims2, s3, level_list, file_mapping):
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
            ds = datasets[(i, j)]
            for name, dim in ds.dimensions.items():
                if (i == 0 and name.startswith(s1)) or (j == 0 and name.startswith(s2)):
                    dst_dims[name] += len(dim)
                else:
                    if dst_dims[name] != len(dim):
                        raise Exception("The blocks don't fit dude!")
    for name, dim in src.dimensions.items():
        dst.createDimension(name, (dst_dims[dim] if not dim.isunlimited() else None))
    if any(level_list):
        dst.createDimension(s3, len(level_list))

    # Copy variables
    dst_vars, time_indices = {}, {}
    for name, variable in src.variables.items():
        dst_vars[name] = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst_vars[name].setncatts(src.variables[name].__dict__)
        time_indices[name] = match_dim_character(name, variable, "time")

    for i in range(len(src.dimensions["time"])):  # loop over t
        var_blocks = {}
        for j in dims2:  # loop over y
            var_slices = {}
            ds = datasets[(0, j)]
            for varname, vardata in ds.variables.items():
                if time_indices[varname] < 0:
                    if i == 0:
                        var_slices[varname] = vardata
                else:
                    var_slices[varname] = vardata.take(indices=i, axis=time_indices[varname])
            for k in dims1:  # loop over x
                ds = datasets[(k, j)]
                for varname, vardata in ds.variables.items():
                    axis_index = match_dim_character(varname, vardata, s1)
                    if axis_index > 0:
                        values = vardata.take(indices=i, axis=time_indices[varname])
                        var_slices[varname] = var_slices[varname].append(values, axis=axis_index)
            for varname, vardata in var_blocks.items():
                axis_index = match_dim_character(varname, vardata, s2)
                if axis_index > 0:
                    var_blocks[varname] = var_blocks[varname].append(var_slices[varname], axis=axis_index)
        for varname, vardata in var_blocks.items():
            time_index = time_indices[varname]
            if time_index >= 0:
                numpy.copyto(dst_vars[varname].take(indices=i, axis=time_index), var_blocks[varname])
            elif i == 0:
                numpy.copyto(dst_vars[varname], var_blocks[varname])


def main():
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("datadir", metavar="DIR", type=str, help="Dales output (run) directory")
    parser.add_argument("--exp", "-e", metavar="N", type=int, default=1, help="Experiment number (default: 001)")
    args = parser.parse_args()
    data_dir = args.datadir
    dims1, dims2, levs, files = find_xsec(data_dir, args.exp, 'x', 'y')
    build_xsec('x', dims1, 'y', dims2, 'z', levs, files)


if __name__ == "__main__":
    main()
