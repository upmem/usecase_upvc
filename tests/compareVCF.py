#!/usr/bin/env python

import sys
import re


def per(val, total):
    return (val * 100.0) / total


stat_limit = 100


def init_stat():
    stat = {}
    for i in range(0, stat_limit):
        stat[i] = [0, 0, 0]
    return stat


def extract_info(info):
    regex = re.search('DEPTH=(.*);COV=(.*);SCORE=(.*)', info)
    depth = int(regex.group(1))
    cov = int(regex.group(2))
    score = int(regex.group(3))
    percentage = 100.0
    if cov != 0:
        percentage = int(per(depth, cov))
    return depth, percentage, score


def update_stat_no_check(stat, v):
    info = v[2]
    depth, percentage, score = extract_info(info)

    if depth >= stat_limit:
        depth = stat_limit - 1
    stat[depth][0] += 1

    if percentage >= stat_limit:
        percentage = stat_limit - 1
    stat[percentage][1] += 1

    if score >= stat_limit:
        score = stat_limit - 1
    stat[score][2] += 1


def update_stat(stat, v):
    if stat is None:
        return
    update_stat_no_check(stat, v)


def update_stat_for_pos_no_check(stat, V):
    for v in V:
        update_stat_no_check(stat, v)


def update_stat_for_pos(stat, V):
    if stat is None:
        return
    update_stat_for_pos_no_check(stat, V)


def update_stat_for_chr(stat, V):
    if stat is None:
        return
    for v in V.values():
        update_stat_for_pos_no_check(stat, v)


def print_stat(tp_stat, fp_stat, nb_tp, nb_fp, nb_cm):
    if tp_stat is None or fp_stat is None:
        return
    if nb_fp == 0:
        nb_fp = -1
    if nb_tp == 0:
        nb_tp = -1
    if nb_cm == 0:
        nb_cm = -1
    tp_acc_depth = 0
    tp_acc_percentage = 0
    tp_acc_score = 0
    fp_acc_depth = 0
    fp_acc_percentage = 0
    fp_acc_score = 0
    print("1%% tp = %.2f%% cm, 1%% cm = %.2f%% fp" %
          ((nb_tp / float(nb_cm)), (nb_cm / float(nb_fp))))
    print("i\t| tp/depth\t\tfp/depth\t| tp/percentage\t\tfp/percentage\t| tp/score\t\tfp/score")
    print("--------"
          "+---------------------------------------"
          "+---------------------------------------"
          "+---------------------------------------"
          )
    for i in range(0, stat_limit):
        tp_acc_depth += tp_stat[i][0]
        tp_acc_percentage += tp_stat[i][1]
        tp_acc_score += tp_stat[i][2]
        fp_acc_depth += fp_stat[i][0]
        fp_acc_percentage += fp_stat[i][1]
        fp_acc_score += fp_stat[i][2]
        print("%d\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)" %
              (i,
               per(tp_acc_depth, nb_tp),
               per(tp_stat[i][0], nb_tp),
               per(fp_acc_depth, nb_fp),
               per(fp_stat[i][0], nb_fp),
               per(tp_acc_percentage, nb_tp),
               per(tp_stat[i][1], nb_tp),
               per(fp_acc_percentage, nb_fp),
               per(fp_stat[i][1], nb_fp),
               100.0 - per(tp_acc_score, nb_tp),
               per(tp_stat[i][2], nb_tp),
               100.0 - per(fp_acc_score, nb_fp),
               per(fp_stat[i][2], nb_fp),
               )
              )


def var_match(v1, v2):
    return v1[0] == v2[0] and v1[1] == v2[1]


def compute_for_pos(V1, V2, tp, fp, tp_stat, fp_stat):
    for v1 in V1:
        found = False
        for v2 in V2:
            if var_match(v1, v2):
                found = True
                tp += 1
                update_stat(tp_stat, v1)
                break
        if not found:
            fp += 1
            update_stat(fp_stat, v1)
    return tp, fp


def compute(V1, V2, tp_stat, fp_stat):
    tp = 0
    fp = 0
    for chr, V1_chr in V1.items():
        if chr not in V2:
            update_stat_for_chr(fp_stat, V1_chr)
            continue
        for pos, V1_pos in V1_chr.items():
            if pos in V2[chr]:
                tp, fp = compute_for_pos(
                    V1_pos, V2[chr][pos], tp, fp, tp_stat, fp_stat)
            else:
                fp += len(V1_pos)
                update_stat_for_pos(fp_stat, V1_pos)
    return tp, fp


def compute_data(V_ref, V_upvc, len_ref, len_upvc):
    enable_stat = False
    tp_stat = None
    fp_stat = None
    if enable_stat:
        tp_stat = init_stat()
        fp_stat = init_stat()

    tp, fp = compute(V_upvc, V_ref, tp_stat, fp_stat)
    cm, fn = compute(V_ref, V_upvc, None, None)

    print_stat(tp_stat, fp_stat, tp, fp, cm)

    if tp != cm:
        print("ERROR while computing TP et CM")
    if tp + fp != len_upvc:
        print("ERROR while computing TP and FP")
    if cm + fn != len_ref:
        print("ERROR while computing CM and FN")

    print("tp:\t%.2f%%\t(%d/%d)" % (per(tp, len_upvc), tp, len_upvc))
    print("fp:\t%.2f%%\t(%d/%d)" % (per(fp, len_upvc), fp, len_upvc))
    print("fn:\t%.2f%%\t(%d/%d)" % (per(fn, len_ref), fn, len_ref))
    print("cm:\t%.2f%%\t(%d/%d)" % (per(cm, len_ref), cm, len_ref))


def filter(info, per_filter, score_filter):
    depth, percentage, score = extract_info(info)
    if score <= score_filter and percentage >= per_filter:
        return True
    return False


def get_data(filename, filter_enable, per_filter, score_filter):
    print("\nreading " + filename)
    SUB = {}
    INS = {}
    DEL = {}
    ff = open(filename)
    line = ff.readline()
    len_sub = 0
    len_ins = 0
    len_del = 0
    while line != "":
        if line[0] == "#":
            line = ff.readline()
            continue
        L = line.split("\t")
        if L[0] == "X":
            chr = 23
        elif L[0] == "Y":
            chr = 24
        else:
            chr = int(re.search('[^0-9]*(.*)', L[0]).group(1))
        pos = int(L[1])
        ref = L[3]
        alt = L[4]
        info = L[7]
        if len(ref) > 1 and len(alt) <= 1 and (not filter_enable or filter(info, per_filter, score_filter)):
            if chr not in DEL:
                DEL[chr] = {}
            if pos not in DEL[chr]:
                DEL[chr][pos] = [[ref, alt, info]]
            else:
                DEL[chr][pos].append([ref, alt, info])
            len_del += 1
        elif len(alt) > 1 and len(ref) <= 1 and (not filter_enable or filter(info, per_filter, score_filter)):
            if chr not in INS:
                INS[chr] = {}
            if pos not in INS[chr]:
                INS[chr][pos] = [[ref, alt, info]]
            else:
                INS[chr][pos].append([ref, alt, info])
            len_ins += 1
        elif len(ref) == 1 and len(alt) == 1 and (not filter_enable or filter(info, per_filter, score_filter)):
            if chr not in SUB:
                SUB[chr] = {}
            if pos not in SUB[chr]:
                SUB[chr][pos] = [[ref, alt, info]]
            else:
                SUB[chr][pos].append([ref, alt, info])
            len_sub += 1
        line = ff.readline()

    ff.close()
    print(len_sub, len_ins, len_del)
    return SUB, INS, DEL, len_sub, len_ins, len_del

ref_file = sys.argv[1]
upvc_file = sys.argv[2]
SUB_ref, INS_ref, DEL_ref, len_ref_sub, len_ref_ins, len_ref_del = get_data(
    ref_file, False, 0, 0)
SUB_upvc, INS_upvc, DEL_upvc, len_upvc_sub, len_upvc_ins, len_upvc_del = get_data(
    upvc_file, False, 0, 0)

print("\nsubstitution")
compute_data(SUB_ref, SUB_upvc, len_ref_sub, len_upvc_sub)

print("\ninsertions")
compute_data(INS_ref, INS_upvc, len_ref_ins, len_upvc_ins)

print("\ndeletion")
compute_data(DEL_ref, DEL_upvc, len_ref_del, len_upvc_del)
