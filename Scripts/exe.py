import os
import subprocess


PATH = "../CPP"


def run(dataset, k = 0, noscheme = False):
    # set the path in which to run the executable
    old_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(PATH)

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
            line = line.split()
            k = line[-3]
            k = k[2:]
        if "Done" in line:
            line = line.split()
            iters = line[-2]
            runtime = line[2]
    return docs, words, nz, k, iters, runtime


def main(fname, outfile):
    # make sure file exists
    if not os.path.isfile(fname):
        print("File \"" + fname + "\" not found.")
        return

    # output the processed data into an output file
    output = open(outfile, "w")
    output.write("DC        WC        NUM_NZ    K         ITERS     RUNTIME\n")
    output.close()

    # read file and process each line
    run_data = []
    runs = open(fname, "r")
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
        out, err = run(dataset, k, noscheme)
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


main("runs.txt", "data.txt")
