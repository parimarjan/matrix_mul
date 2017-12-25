import subprocess as sp
import time
import argparse
from collections import defaultdict
import re

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

def run_matrix_test(num_els):
    cmd_template = './test -n {n}'
    cmd = cmd_template.format(n = num_els)
    cmd = cmd.split()
    timings = defaultdict(list)
    avg_timings = defaultdict(float)

    for i in range(NUM_TRIES):
        print('iteration: ', i)
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
args = parser.parse_args()

power_timings = []
for power in range(args.start,args.end,1):
    num_els = 2**power
    power_timings.append(run_matrix_test(num_els))

# plotting time!
for k in power_timings[0]:
    # each of the keys, like naive etc.
    timings = [p[k] for p in power_timings]

for t in timings:
    print(t)


