#!/usr/bin/env python

import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

if len(sys.argv) != 2:
    sys.stderr.write("Usage: cat windows.bed | ./calc_window_mappability.py mappability.bed\n")
    sys.exit(1)

window_file = sys.stdin
values_file = open(sys.argv[1])
last_val    = None

def read_window():
    line = window_file.readline()
    if len(line) == 0:
        return None
    i1 = line.find('\t')
    i2 = line.find('\t', i1+1)
    i3 = line.find('\t', i2+1);
    if (i3 == -1):
        i3 = len(line) - 1;
    return (line[:i1], int(line[i1+1:i2]), int(line[i2+1:i3]), line[i3:-1])

def read_value():
    line = values_file.readline()
    if len(line) == 0:
        return None
    i1 = line.find('\t')
    i2 = line.find('\t', i1+1)
    i3 = line.find('\t', i2+1)
    return (line[:i1], int(line[i1+1:i2]), int(line[i2+1:i3]), float(line[i3+1:-1]))

while True:
    window = read_window()
    if window is None:
        break

    tot_val = 0.0
    tot_bp  = 0
    first = True

    while True:
        if first and last_val is not None and last_val[0] == window[0]:
            value = last_val
        else:
            value = read_value()
            last_val = value
        first = False

        if value is None:
            break
        while value[0] != window[0]:
            value = read_value()
            last_val = value
        if value[2] <= window[1]:
            continue
        if value[1] >= window[2]:
            break

        n_bp = min(value[2], window[2]) - max(value[1], window[1])
        tot_val += n_bp * value[3]
        tot_bp  += n_bp

        if value[2] >= window[2]:
            break

    print "%s\t%d\t%d%s\t%.3f" % (window[0], window[1], window[2], window[3],
                                  tot_val/tot_bp if tot_bp > 0 else 0.0)

window_file.close()
values_file.close()

