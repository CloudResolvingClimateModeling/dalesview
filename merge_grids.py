#!/usr/bin/env python

import argparse
import os
import re


def find_xsec(data_dir, expnr, s1, s2, s3):
    regex = re.compile("^cross" + s1 + s2 + ".(\d+)." + s1 + "(\d+)" + s2 + "(\d+)." + str(expnr).zfill(3) + ".nc$")
    for filepath in os.listdir(data_dir):
        if re.match(regex, os.path.basename(filepath)):
            result = re.search(regex, os.path.basename(filepath))
            tuple = (int(result.group(2)), int(result.group(3)), float(result.group(1)))
            print "Found file:" + filepath + " matches are " + str(tuple)



def main():
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("datadir", metavar="DIR", type=str, help="Dales output (run) directory")
    parser.add_argument("--exp", "-e", metavar="N", type=int, default=1, help="Experiment number (default: 001)")
    args = parser.parse_args()
    data_dir = args.datadir
    find_xsec(data_dir, args.exp, 'x', 'y', 'z')


if __name__ == "__main__":
    main()
