#!/usr/bin/env python

import sys
import re

stat_limit = 100
def init_stat():
    stat = {}
    for i in range (0, stat_limit):
        stat[i] = [0, 0]
    return stat

def per(val, total):
    return (val * 100) / total
def per_str(val, total):
    return str(per(val, total)) + "%"

def extract_info(info):
    regex = re.search('DEPTH=(.*);COV=(.*)', info)
    depth = int(regex.group(1), 10)
    cov = int(regex.group(2), 10)
    percentage = 0
    if cov != 0:
        percentage = per(depth, cov)
    return depth, cov, percentage

def update_stat(stat, info):
    depth, cov, percentage = extract_info(info)
    if depth >= stat_limit:
        stat[stat_limit - 1][0] += 1
    else:
        stat[depth][0] += 1
    if percentage >= stat_limit:
        stat[stat_limit - 1][1] += 1
    elif percentage < 0:
        print(info, depth, cov, percentage)
    else:
        stat[percentage][1] += 1

def print_stat(tp_stat, fp_stat, total):
    tp_acc_depth = 0
    tp_acc_percentage = 0
    fp_acc_depth = 0
    fp_acc_percentage = 0
    for i in range (0, stat_limit):
        tp_acc_depth += tp_stat[i][0]
        tp_acc_percentage += tp_stat[i][1]
        fp_acc_depth += fp_stat[i][0]
        fp_acc_percentage += fp_stat[i][1]
        print(str(i) + "\t"
              + "tp_d:\t" + str(tp_stat[i][0]) + "\t" + per_str(tp_stat[i][0] , total) + "\t"
              + "tp_p:\t" + str(tp_stat[i][1]) + "\t" + per_str(tp_stat[i][1] , total) + "\t"
              + "fp_d:\t" + str(fp_stat[i][0]) + "\t" + per_str(fp_stat[i][0] , total) + "\t"
              + "fp_p:\t" + str(fp_stat[i][1]) + "\t" + per_str(fp_stat[i][1] , total) + "\t"
              + "tp_d: " + per_str(tp_acc_depth, total) + "\t"
              + "tp_p: " + per_str(tp_acc_percentage, total) + "\t"
              + "fp_d: " + per_str(fp_acc_depth, total) + "\t"
              + "fp_p: " + per_str(fp_acc_percentage, total) + "\t"
              )

delta = 10

def var_match(V1, V2):
    return V1[0] == V2[0] and V1[1] == V2[1]

def filter(info, f_ok, f_next, d_ok, d_next, p_ok, p_next):
    depth, cov, percentage = extract_info(info)
    if depth >= d_ok and percentage >= p_ok:
        f_ok += 1
    elif depth >= d_next and percentage >= p_next:
        f_next += 1
    return f_ok, f_next


def computeTPFP(V_ref, V_upvc, filename_prefix, d_ok, d_next, p_ok, p_next):
    file_tp = open(filename_prefix + ".tp", "w")
    file_fp = open(filename_prefix + ".fp", "w")
    tp = 0
    fp = 0
    tp_f_ok = 0
    tp_f_next = 0
    fp_f_ok = 0
    fp_f_next = 0
    tp_stat = init_stat()
    fp_stat = init_stat()
    for i in V_upvc :
        found = False
        for d in range(-delta, delta):
            if i + d in V_ref and var_match(V_upvc[i], V_ref[i + d]):
                tp_f_ok, tp_f_next = filter(V_upvc[i][2], tp_f_ok, tp_f_next, d_ok, d_next, p_ok, p_next)
                tp += 1
                file_tp.write(str(V_upvc[i]) + "\n")
                update_stat(tp_stat, V_upvc[i][2])
                found = True
                break
        if not found:
            fp_f_ok, fp_f_next = filter(V_upvc[i][2], fp_f_ok, fp_f_next, d_ok, d_next, p_ok, p_next)
            fp += 1
            file_fp.write(str(V_upvc[i]) + "\n")
            update_stat(fp_stat, V_upvc[i][2])
    print_stat(tp_stat, fp_stat, len(V_upvc))
    print("After filtering: ", tp_f_ok, tp_f_next, tp - tp_f_ok - tp_f_next, fp_f_ok, fp_f_next)
    file_tp.close()
    file_fp.close()
    return tp, fp

def computeFN(V_ref, V_upvc):
    fn = 0

    for i in V_ref:
        found = False
        for d in range (-delta, delta):
            if  i + d in V_upvc and var_match(V_upvc[i + d], V_ref[i]):
                found = True
                break;
        if not found:
            fn += 1

    return fn

def computeData (V_ref, V_upvc, filename_prefix, d_ok, d_next, p_ok, p_next):
    tp , fp = computeTPFP(V_ref, V_upvc, filename_prefix, d_ok, d_next, p_ok, p_next)
    fn = computeFN(V_ref, V_upvc)

    print ("tp:\t" + per_str(tp, len(V_upvc)) + "\t(" + str(tp) + " / " + str(len(V_upvc)) + ")")
    print ("fp:\t" + per_str(fp, len(V_upvc)) + "\t(" + str(fp) + " / " + str(len(V_upvc)) + ")")
    print ("fn:\t" + per_str(fn, len(V_ref))  + "\t(" + str(fn) + " / " + str(len(V_ref)) + ")")
    print ("cm:\t" + per_str(tp, len(V_ref))  + "\t(" + str(tp) + " / " + str(len(V_ref)) + ")")

def getData(filename):
    print("\nreading " + filename)
    SUB = {}
    INS = {}
    DEL = {}
    ff = open(filename)
    line = ff.readline()
    while line != "":
        if line[0] == "#":
            line = ff.readline()
            continue
        L = line.split("\t")
        pos = int(L[1])
        ref = L[3]
        alt = L[4]
        info = L[7]
        if len(ref) > 1 and len(alt) == 1:
            DEL[pos] = [ref,alt,info]
        elif len(alt) > 1 and len(ref) == 1:
            INS[pos] = [ref,alt,info]
        elif len(ref) == 1 and len(alt) == 1:
            SUB[pos] = [ref,alt,info]
        else:
            print("!!" + line[0:-1])
        line = ff.readline()

    ff.close()
    print (len(SUB), len(INS), len(DEL))
    return SUB, INS, DEL

ref_file = sys.argv[1]
upvc_file = sys.argv[2]
SUB_ref, INS_ref, DEL_ref = getData(ref_file)
SUB_upvc, INS_upvc, DEL_upvc = getData(upvc_file)

print ("\nsubstitution")
computeData(SUB_ref, SUB_upvc, upvc_file, 3, 2, 30, 15)

print ("\ninsertions")
computeData(INS_ref, INS_upvc, upvc_file, 2, 1, 17, 2)

print ("\ndeletion")
computeData(DEL_ref, DEL_upvc, upvc_file, 2, 1, 17, 2)


