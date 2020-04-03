#!/usr/bin/env python

import sys
import re

stat_limit = 50
def init_stat():
    stat = {}
    for i in range (0, stat_limit):
        stat[i] = [0, 0, 0]
    return stat

def per(val, total):
    return (val * 100.0) / total
def per_str(val, total):
    return str(per(val, total))[0:5] + "%"

def extract_info(info):
    regex = re.search('DEPTH=(.*);COV=(.*);SCORE=(.*)', info)
    depth = int(regex.group(1), 10)
    cov = int(regex.group(2), 10)
    score = int(regex.group(3), 10)
    average_score = score / depth
    percentage = 0.0
    if cov != 0:
        percentage = per(depth, cov)
    return depth, cov, percentage, average_score

def update_stat(stat, depth, percentage, score):
    if depth >= stat_limit:
        stat[stat_limit - 1][0] += 1
    else:
        stat[depth][0] += 1
    if percentage >= stat_limit:
        stat[stat_limit - 1][1] += 1
    else:
        stat[int(percentage)][1] += 1
    if score >= stat_limit:
        stat[stat_limit - 1][2] += 1
    else:
        stat[score][2] += 1

def print_stat(tp_stat, fp_stat, tp_total, upvc_total):
    tp_acc_depth = 0
    tp_acc_percentage = 0
    tp_acc_score = 0
    fp_acc_depth = 0
    fp_acc_percentage = 0
    fp_acc_score = 0
    for i in range (0, stat_limit):
        tp_acc_depth += tp_stat[i][0]
        tp_acc_percentage += tp_stat[i][1]
        tp_acc_score += tp_stat[i][2]
        fp_acc_depth += fp_stat[i][0]
        fp_acc_percentage += fp_stat[i][1]
        fp_acc_score += fp_stat[i][2]
        print(str(i) + "\t"
              + "tp_d: " + per_str(tp_stat[i][0] , tp_total) + "\t"
              + "tp_p: " + per_str(tp_stat[i][1] , tp_total) + "\t"
              + "tp_s: " + per_str(tp_stat[i][2] , tp_total) + "\t"
              + "fp_d: " + per_str(fp_stat[i][0] , upvc_total) + "\t"
              + "fp_p: "+ per_str(fp_stat[i][1] , upvc_total) + "\t"
              + "fp_s: " + per_str(fp_stat[i][2] , upvc_total) + "\t"
              + "tp_d: " + per_str(tp_acc_depth, tp_total) + "\t"
              + "tp_p: " + per_str(tp_acc_percentage, tp_total) + "\t"
              + "tp_s: " + per_str(tp_acc_score, tp_total) + "\t"
              + "fp_d: " + per_str(fp_acc_depth, upvc_total) + "\t"
              + "fp_p: " + per_str(fp_acc_percentage, upvc_total) + "\t"
              + "fp_s: " + per_str(fp_acc_score, upvc_total) + "\t"
              )

def var_match(V1, V2):
    return V1[0] == V2[0] and V1[1] == V2[1]

def filter(res, depth, percentage, score, limit_d1, limit_d2, limit_p, limit_s):
    if depth >= limit_d1:
        return res + 1
    if depth >= limit_d2 and percentage >= limit_p and score <= limit_s:
        return res + 1
    return res

def computeTPFP(V_ref, V_upvc, isDel, limit_d1, limit_d2, limit_p, limit_s):
    tp = 0
    fp = 0
    tp_f = 0
    fp_f = 0
    tp_stat = init_stat()
    fp_stat = init_stat()
    for i in V_upvc :
        depth, cov, percentage, score = extract_info(V_upvc[i][2])
        d = 0
        if isDel:
            d = -len(V_upvc[i][0])
        if i + d in V_ref and var_match(V_upvc[i], V_ref[i + d]):
            tp += 1
            tp_f = filter(tp_f, depth, percentage, score, limit_d1, limit_d2, limit_p, limit_s)
            update_stat(tp_stat, depth, percentage, score)
        else:
            fp += 1
            fp_f = filter(fp_f, depth, percentage, score, limit_d1, limit_d2, limit_p, limit_s)
            update_stat(fp_stat, depth, percentage, score)
    print_stat(tp_stat, fp_stat, tp, len(V_upvc))
    print("Filtered:", tp_f, fp_f)
    return tp, fp

def computeFN(V_ref, V_upvc, isDel):
    fn = 0

    for i in V_ref:
        d = 0
        if isDel:
            d = len(V_ref[i][0])
        if not i + d in V_upvc or not var_match(V_upvc[i + d], V_ref[i]):
            fn += 1

    return fn

def computeData (V_ref, V_upvc, isDel, limit_d1, limit_d2, limit_p, limit_s):
    tp , fp = computeTPFP(V_ref, V_upvc, isDel, limit_d1, limit_d2, limit_p, limit_s)
    fn = computeFN(V_ref, V_upvc, isDel)

    if tp + fp != len(V_upvc):
        print("ERROR while computing TP and FP")
    if tp + fn != len(V_ref):
        print("ERROR while computing FN")

    print ("tp:\t" + per_str(tp, len(V_upvc)) + "\t(" + str(tp) + " / " + str(len(V_upvc)) + ")")
    print ("fp:\t" + per_str(fp, len(V_upvc)) + "\t(" + str(fp) + " / " + str(len(V_upvc)) + ")")
    print ("fn:\t" + per_str(fn, len(V_ref))  + "\t(" + str(fn) + " / " + str(len(V_ref))  + ")")
    print ("cm:\t" + per_str(tp, len(V_ref))  + "\t(" + str(tp) + " / " + str(len(V_ref))  + ")")

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
        if len(ref) > 1 and len(alt) <= 1:
            DEL[pos] = [ref, alt, info]
        elif len(alt) > 1 and len(ref) <= 1:
            INS[pos] = [ref, alt, info]
        elif len(ref) == 1 and len(alt) == 1:
            SUB[pos] = [ref, alt, info]
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
computeData(SUB_ref, SUB_upvc, False, 2, 1, 2, 29)

print ("\ninsertions")
computeData(INS_ref, INS_upvc, False, 4, 2, 15, 19)

print ("\ndeletion")
computeData(DEL_ref, DEL_upvc, True, 4, 2, 15, 19)
