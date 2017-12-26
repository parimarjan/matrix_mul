import subprocess as sp
import time
import argparse
from collections import defaultdict
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json

NUM_TRIES = 5

def get_time_from_string(string):
    '''
    '''
    # matches scientific notation and stuff.
    numbers = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",
            string)
    # all the things we are printing so far has just 1 num.
    if len(numbers) >= 1:
        return float(numbers[0])
    else:
        return None

def run_matrix_test(num_els, cmd_template):
    cmd = cmd_template.format(n = num_els)
    cmd = cmd.split()
    timings = defaultdict(list)
    avg_timings = defaultdict(float)

    for i in range(NUM_TRIES):
        process = sp.Popen(cmd, stdout=sp.PIPE)
        # process.wait()
        stdout = process.communicate()[0]
        stdout = stdout.split('\n')
        for s in stdout:
            if "time" in s:
                # one of the timing tests.
                name = s.split()[0]
                timings[name].append(get_time_from_string(s))

    for k,v in timings.iteritems(): 
        avg_timings[k] = sum(v) / len(v)

    return avg_timings

parser = argparse.ArgumentParser()
parser.add_argument("-start", "-s", type=int, required=False,
                    default=8, help="start index of powers")
parser.add_argument("-end", "-e", type=int, required=False,
                    default=9, help="end index of powers")
parser.add_argument("-file", "-f", type=str, required=False,
                    default="fig.png", help="name of output figure")

args = parser.parse_args()

# power_timings = []

all_timings = defaultdict(list)

c_template = './test -n {n}'
js_template = 'node test.js -n {n}'

for power in range(args.start,args.end,1):
    num_els = 2**power
    avg_timings = run_matrix_test(num_els, c_template)
    for k,v in avg_timings.iteritems():
        all_timings[k].append(v)

    avg_timings2 = run_matrix_test(num_els, js_template)
    for k, v in avg_timings2.iteritems():
        all_timings['javascript'].append(float(v)/1000)


x_vals = [2**p for p in range(args.start, args.end, 1)]
x_label = 'dimension of multiplied square matrices'
y_label = 'time in seconds'
colors = ["r-", "b--", "g-", "c--"]

i = 0
for k,v in all_timings.iteritems():
    # dump a json with this k
    fname = k + '.json'
    data = []
    for j,val in enumerate(v): 
        data.append([x_vals[j], val])
    with open(fname, 'w') as outfile:
        json.dump(data, outfile)

    plt.plot(x_vals, v, colors[i], label=k)
    i += 1

plt.xlabel(x_label, size=15)
plt.ylabel(y_label, size=15)
plt.legend(loc='best')
# plt.show()
plt.savefig(args.file)


