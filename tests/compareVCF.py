#!/usr/bin/env python3

import sys
import re


def per(val, total):
    return (val * 100.0) / total


stat_limit = 100


def init_stat():
    return [[0, 0, 0] for _ in range(stat_limit)]


extract_info_re = re.compile("DEPTH=(.*);COV=(.*);SCORE=(.*)")


def extract_info(info):
    regex = re.search(extract_info_re, info)
    depth = int(regex.group(1))
    cov = int(regex.group(2))
    score = int(regex.group(3))
    percentage = 100.0
    if cov != 0:
        percentage = int(per(depth, cov))
    return depth, percentage, score


def update_stat_no_check(stat, info):
    depth, percentage, score = extract_info(info)

    stat[depth if depth < stat_limit else stat_limit - 1][0] += 1
    stat[percentage if percentage < stat_limit else stat_limit - 1][1] += 1
    stat[score if score < stat_limit else stat_limit - 1][2] += 1


def update_stat(stat, info):
    if stat is None:
        return
    update_stat_no_check(stat, info)


def update_stat_for_pos(stat, infos):
    if stat is None:
        return
    for info in infos.values():
        update_stat_no_check(stat, info)


def print_stat(tp_stat, fp_stat, nb_tp, nb_fp, nb_cm):
    if tp_stat is None or fp_stat is None:
        return

    tp_acc_depth = 0
    tp_acc_percentage = 0
    tp_acc_score = nb_tp
    fp_acc_depth = 0
    fp_acc_percentage = 0
    fp_acc_score = nb_fp

    if nb_fp == 0:
        nb_fp = -1
    if nb_tp == 0:
        nb_tp = -1
    if nb_cm == 0:
        nb_cm = -1

    print("1%% tp = %.2f%% cm, 1%% cm = %.2f%% fp" %
          ((nb_tp / float(nb_cm)), (nb_cm / float(nb_fp))))
    print("i\t| tp/depth\t\tfp/depth\t| tp/percentage\t\tfp/percentage\t| tp/score\t\tfp/score")
    print("--------"
          "+---------------------------------------"
          "+---------------------------------------"
          "+---------------------------------------"
          )
    for i in range(stat_limit):
        tp_depth = tp_stat[i][0]
        tp_per = tp_stat[i][1]
        tp_score = tp_stat[i][2]
        fp_depth = fp_stat[i][0]
        fp_per = fp_stat[i][1]
        fp_score = fp_stat[i][2]
        tp_acc_depth += tp_depth
        tp_acc_percentage += tp_per
        tp_acc_score -= tp_score
        fp_acc_depth += fp_depth
        fp_acc_percentage += fp_per
        fp_acc_score -= fp_score
        print("%d\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)\t| %.2f%% (%.2f%%)   \t%.2f%% (%.2f%%)" %
              (i,
               per(tp_acc_depth, nb_tp), per(tp_depth, nb_tp),
               per(fp_acc_depth, nb_fp), per(fp_depth, nb_fp),
               per(tp_acc_percentage, nb_tp), per(tp_per, nb_tp),
               per(fp_acc_percentage, nb_fp), per(fp_per, nb_fp),
               per(tp_acc_score, nb_tp), per(tp_score, nb_tp),
               per(fp_acc_score, nb_fp), per(fp_score, nb_fp),
               )
              )


def compute_for_pos(v1_infos, v2_infos, tp, fp, tp_stat, fp_stat):
    for (ref, alt), v1_info in v1_infos.items():
        if (ref, alt) in v2_infos:
            tp += 1
            update_stat(tp_stat, v1_info)
        else:
            fp += 1
            update_stat(fp_stat, v1_info)
    return tp, fp


def compute(V1, V2, tp_stat, fp_stat):
    tp = 0
    fp = 0
    for (chr, pos), v1_infos in V1.items():
        if (chr, pos) in V2:
            tp, fp = compute_for_pos(
                v1_infos, V2[(chr, pos)], tp, fp, tp_stat, fp_stat)
        else:
            fp += len(v1_infos)
            update_stat_for_pos(fp_stat, v1_infos)
    return tp, fp


def compute_data(V_ref, V_upvc, len_ref, len_upvc):
    enable_stat = False
    tp_stat = init_stat() if enable_stat else None
    fp_stat = init_stat() if enable_stat else None

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
    return score <= score_filter and percentage >= per_filter


def get_data(filename, filter_enable, per_filter, score_filter):
    print("\nreading " + filename)
    SUB = {}
    INS = {}
    DEL = {}
    len_sub = 0
    len_ins = 0
    len_del = 0

    for line in open(filename, "r").readlines():
        if len(line) < 1 or line[0] == "#":
            continue

        chr_str, pos_str, _, ref, alt, _, _, info = line.split("\t")

        if chr_str == "X":
            chr = 23
        elif chr_str == "Y":
            chr = 24
        else:
            chr = int(chr_str.strip("chr"))

            pos = int(pos_str)

            if len(ref) > 1 and len(alt) <= 1 \
               and (not filter_enable or filter(info, per_filter, score_filter)):
                DEL.setdefault((chr, pos), {})[(ref, alt)] = info
                len_del += 1
            elif len(alt) > 1 and len(ref) <= 1 \
                 and (not filter_enable or filter(info, per_filter, score_filter)):
                INS.setdefault((chr, pos), {})[(ref, alt)] = info
                len_ins += 1
            elif len(ref) == 1 and len(alt) == 1 \
                 and (not filter_enable or filter(info, per_filter, score_filter)):
                SUB.setdefault((chr, pos), {})[(ref, alt)] = info
                len_sub += 1

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
