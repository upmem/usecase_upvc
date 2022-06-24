#!/usr/bin/env python3


import sys


nb_dpu = 0
last_dpu_id = -1
dict = {}

for line in open(sys.argv[1], "r").readlines():
    line_elems = line.split()
    log = line_elems[0]
    if log != "LOG":
        continue
    dpu_id = int(line_elems[1].split("=")[1])
    if dpu_id != last_dpu_id:
        nb_dpu += 1
        last_dpu_id = dpu_id

    key, value = line_elems[2].split("=")
    dict[key] = dict.setdefault(key, 0) + int(value)

print("nb_dpu:", nb_dpu)
print(dict)

for key, value in dict.items():
    dict[key] = value / nb_dpu

print(dict)
