#!/usr/bin/env python

import argparse
import itertools
import multiprocessing
import os
import re
import netCDF4
import numpy
import logging

log = logging.getLogger(__name__)

chunk_size = 10
n_digits = 6


def merge_files(data_dir, regex, output, expnr):
    file_info = find_files(data_dir, regex)
    for exp, (keyvals, filemap) in file_info.items():
        if 0 < expnr != int(exp):
            continue
        if any(filemap):
            dims1 = keyvals[0]
            dims2 = keyvals[1] if len(keyvals) > 1 else [1]
            levs = keyvals[2] if len(keyvals) > 2 else []
            glue_grids(dims1, dims2, levs, filemap, '.'.join([output, exp, "nc"]))


def find_files(data_dir, regex_str):
    regex = re.compile(regex_str)
    axes = ("x", "y", "lev")
    file_mappings = {}
    for filepath in os.listdir(data_dir):
        if re.match(regex, os.path.basename(filepath)):
            result = re.search(regex, os.path.basename(filepath)).groupdict()
            key = []
            for axis in axes:
                if result.get(axis, None) is not None:
                    key.append(result[axis])
            key = tuple(key)
            exp = result.get("exp", None)
            if exp is None:
                log.warning("Could not detect experiment number for file %s...setting it to 001" % filepath)
                exp = "001"
            if exp in file_mappings:
                file_mappings[exp][key] = os.path.join(data_dir, filepath)
            else:
                file_mappings[exp] = {key: os.path.join(data_dir, filepath)}
    key_values = {}
    for exp, file_mapping in file_mappings.items():
        found_keys = file_mapping.keys()
        if not any(found_keys):
            key_values[exp] = []
        num_keys = len(found_keys[0])
        key_vals = [sorted(list(set([k[i] for k in found_keys]))) for i in range(num_keys)]
        key_values[exp] = key_vals
        for key in itertools.product(*key_vals):
            if key not in found_keys:
                missing_file = regex_str[1:-1]
                for i in range(len(key)):
                    missing_file = missing_file.replace("(?P<" + axes[i] + ">\d+)", key[i])
                missing_file = missing_file.replace("(?P<exp>\d+)", exp)
                log.error("File not found: %s" % missing_file)
                key_values[exp] = []
    result = {}
    for exp, file_mapping in file_mappings.items():
        key_vals = key_values.get(exp, [])
        if any(key_vals):
            new_mapping = {}
            for key, value in file_mapping.items():
                new_key = tuple([int(s) for s in key])
                new_mapping[new_key] = value
            new_key_vals = []
            for values_list in key_vals:
                new_key_vals.append(sorted([int(s) for s in values_list]))
            result[exp] = (new_key_vals, new_mapping)
    return result


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


def glue_grids(dims1, dims2, level_list, file_mapping, output_file):
    global chunk_size, n_digits

    dst = netCDF4.Dataset(output_file, 'w')
    output = os.path.basename(output_file)

    # Copy attributes
    datasets, src = {}, None
    for k in sorted(file_mapping.keys()):
        datasets[k] = netCDF4.Dataset(file_mapping[k], 'r')
        if src is None and len(datasets[k].dimensions.get("time", [])) > 0:
            src = datasets[k]
    if src is None:
        src = datasets[sorted(datasets.keys())[0]]
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

    num_steps = len(src.dimensions["time"])
    num_steps = min(num_steps, len(src.dimensions["time"]))
    dt = min(chunk_size, len(src.dimensions["time"]))

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
        dst_vars[name] = dst.createVariable(name, variable.datatype, dims, zlib=True, least_significant_digit=n_digits)
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
        log.info("processing time step %d of %d for output %s..." % (i, num_steps, output))
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
                    dst_vars[varname][istart:iend, ...] = numpy.swapaxes(time_slices[varname], 0, 1)[istart:iend, ...]
                else:
                    dst_vars[varname][istart:iend, ...] = time_slices[varname][istart:iend, ...]
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
            j = i - 1 if 0 <= z_axis < i else i
            slices.append(slice(0, src.shape[j], 1))
    dest[tuple(slices)] = src[...]


def main():
    global chunk_size, n_digits
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("--dir", metavar="DIR", type=str, default=".", help="Dales output (run) directory")
    parser.add_argument("--odir", metavar="DIR", type=str, default=None,
                        help="Script output directory, by default the run directory")
    parser.add_argument("--exp", "-e", metavar="N", type=int, default=-1, help="Experiment number (default: all)")
    parser.add_argument("--np", "-j", metavar="N", type=int, default=4, help="Number of parallel processes")
    parser.add_argument("--chunksize", metavar="N", type=int, default=10, help="Nr of time slices in memory")
    parser.add_argument("--digits", metavar="N", type=int, default=6,
                        help="Nr of siginificant digits in compressed output")

    args = parser.parse_args()

    data_dir = args.dir
    exp = args.exp
    chunk_size = max(args.chunksize, 1)
    n_digits = max(args.digits, 1)
    n_procs = args.np

    if args.odir is None:
        outdir = data_dir
    else:
        outdir = args.odir
        os.makedirs(outdir)

    dalesfiles = {"crossxy2d": "^crossxy.x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "crossyz2d": "^crossyz.x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "crossxz2d": "^crossxz.x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "crossxy3d": "^crossxy.(?P<lev>\d+).x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "crossyz3d": "^crossyz.(?P<lev>\d+).x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "crossxz3d": "^crossxz.(?P<lev>\d+).x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "surf_xy": "^surf_xy.x(?P<x>\d+)y(?P<y>\d+).(?P<exp>\d+).nc$",
                  "fielddump": "^fielddump.(?P<x>\d+).(?P<y>\d+).(?P<exp>\d+).nc$"}

    def process(ofile):
        merge_files(data_dir, dalesfiles[ofile], os.path.join(outdir, ofile), exp)

    pool = multiprocessing.Pool(processes=n_procs)
    pool.map(process, dalesfiles.keys())


logging.basicConfig(level=logging.DEBUG)

if __name__ == "__main__":
    main()
