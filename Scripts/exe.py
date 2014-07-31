#!/usr/bin/env python3

import os
import subprocess
import sys


EXE_PATH = "../CPP"
IN_FILE = "batch"
OUT_FILE = "output"


def run(path, dataset, k = 0, noscheme = False):
    # set the path in which to run the executable
    old_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(path)

    # set up the value of k
    if k == 0:
        k_cmd = "--autok"
    else:
        k_cmd = "-k " + str(k)

    # set up the command
    cmd = "./spkmeans -d ../TestData/" + dataset + " " + k_cmd + " --noresults --openmp"
    if noscheme:
        cmd += " --noscheme"

    # run the k-means program
    print("Running: \"" + cmd + "\"")
    cmd = cmd.split()
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, \
                              stderr = subprocess.STDOUT)
    out, err = p.communicate()
    if out:
        out = out.decode(encoding = "UTF-8")
    if err:
        err = err.decode(encoding = "UTF-8")

    # change back the path
    os.chdir(old_path)

    return out, err


def interpret(out):
    docs = ""
    words = ""
    nz = ""
    k = ""
    iters = ""
    runtime = ""
    lines = out.split("\n")
    for line in lines:
        if "DATA" in line:
            line = line.split()
            docs = line[1]
            words = line[3]
            nz = line[5]
            nz = nz[1:]
        if "k=" in line:
            if "threads" in line:
                indx = -4
            else:
                indx = -3
            line = line.split()
            k = line[indx]
            k = k[2:]
        if "Done" in line:
            line = line.split()
            iters = line[-2]
            runtime = line[2]
    return docs, words, nz, k, iters, runtime


def main(infile, outfile, path):
    # make sure file exists
    if not os.path.isfile(infile):
        print("File \"" + infile + "\" not found.")
        return

    # output the processed data into an output file
    output = open(outfile, "w")
    output.write("DC        WC        NUM_NZ    K         ITERS     RUNTIME\n")
    output.close()

    # read file and process each line
    run_data = []
    runs = open(infile, "r")
    for line in runs:
        line = line.strip()
        line = line.split()
        if len(line) == 0:
            continue
        dataset = line[0]
        k = 0
        noscheme = False
        if len(line) > 1 and line[1].isdigit():
            k = int(line[1])
        if len(line) > 2 and line[2].lower() == "true":
            noscheme = True
        out, err = run(path, dataset, k, noscheme)
        # if we got anything, format it and append it to the output file
        if out and not err and "error" not in out.lower():
            info = interpret(out)
            output = open(outfile, "a")
            split = list(map(lambda x : x.ljust(10), info[:-1]))
            info = ''.join(split) + info[-1]
            output.write(info + "\n")
            output.close()
            run_data.append(info)
        else:
            print("    ERROR (possibly invalid dataset). Skipping.")
    runs.close()

    # print out the results for all runs
    print("DC        WC        NUM_NZ    K         ITERS     RUNTIME")
    for info in run_data:
        print(info)


if __name__ == "__main__":
    # parse arguments and try to run
    if len(sys.argv) > 1:
        infile = sys.argv[1]
    else:
        infile = IN_FILE
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        outfile = OUT_FILE
    if len(sys.argv) > 3:
        path = sys.argv[3]
    else:
        path = EXE_PATH
    main(infile, outfile, path)
