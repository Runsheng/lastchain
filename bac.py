#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/24 14:17
# @Author  : Runsheng     
# @File    : bac.py

"""
bac,
- to accommodate the specific-designed bac-library, generate a .tab file for the using in sspace
- get bac.tab (from pair-seq) and bac_short.tab (from single-end) for the use in the sspace

start with bam file, so the alignment wrapper is not in this file

WARNING!!!!!
ONLY used for the bac-end data generated by the pBACode, not generalized!
See reference [Wei, 2016] as reference
"pBACode: a random-barcode-based high-throughput approach for BAC paired-end sequencing and physical clone mapping"
"""
from operator import itemgetter
from collections import OrderedDict
import pysam
from Utils import reverse_complement

from longscaf import chain_anchor


def get_tagpair(filename):
    """
    store as a dict
    {"forward":[reverse, number]}
    """

    tagdict={}
    with open(filename, "r") as f:
        for line in f.readlines():
            k=line.split(":")[0]
            v=(line.split(":")[1].strip()).split("\t")
            tagdict[k]=v
    return tagdict


def read_position_store(samfile):
    position_dic = {}

    with pysam.AlignmentFile(samfile, "rb") as samfile:
        for read in samfile.fetch():
            chro = samfile.getrname(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            seq = read.query_sequence
            read_1_or_2 = "1" if read.is_read1 else "2"
            read_F_or_R = "R" if read.is_reverse else "F"
            try:
                position_dic[read.query_name].append((chro, start, end, read_1_or_2, read_F_or_R))
            except KeyError:
                position_dic[read.query_name] = [(chro, start, end, read_1_or_2, read_F_or_R)]

    return position_dic


def read_position_merge(pos_dict, length_cutoff=30):
    """
    merge the forward and reverse reads into one object
    also, get the Barcode
    the new dict as {readname: [barcode, chro, start,end, direction]}
    """
    new_dict = {}
    for k, v in pos_dict.items():
        if len(v) == 2:  # the reads with ambi mapping was ignored
            barcode = k.split("BARCODE")[1]

            chro, start, end, read_1_or_2, read_F_or_R = v[0]
            chro_a, start_a, end_a, read_1_or_2_a, read_F_or_R_a = v[1]
            try:
                if abs(start - end) > length_cutoff and abs(start_a - end_a) > length_cutoff:
                    pos = (start, end, start_a, end_a)
                    if chro == chro_a and read_1_or_2 != read_1_or_2_a and read_F_or_R != read_F_or_R_a:
                        chro_new = chro
                        start_new = min(pos)
                        end_new = max(pos)
                        direction = read_F_or_R if read_1_or_2 == "1" else read_F_or_R_a

                        new_dict[k] = [barcode, chro_new, start_new, end_new, direction]

            except TypeError:  # the break pairs will arise TypeError in (start-end), need to be ignored
                pass

    return new_dict


def read_position_barcode(pos_m):
    """
    pos_l_m, pos_r_m
    change the key to barcode
    the read with same barcode were merged
    data str {barcode: [readname, chro, start,end,direction]}
    """
    barcode_dict = {}

    for k, v in pos_m.items():
        barcode, chro_new, start_new, end_new, direction = v
        try:
            barcode_dict[barcode].append((k, chro_new, start_new, end_new, direction))
        except KeyError:
            barcode_dict[barcode] = []
            barcode_dict[barcode].append((k, chro_new, start_new, end_new, direction))

    return barcode_dict


def glue_FR(tagdict, l_barcode, r_barcode):
    """
    use the tagdict to link the left end and right end

    data str: {barcode_f:
                [(readname1, chro, start,end,direction),
                (readname2, chro, start,end,direction)...],
                [(readname1, chro, start,end,direction)
                (readname2, chro, start,end,direction)...]
              }
    """

    glue_dic = OrderedDict()

    for k, v in tagdict.items():
        barcode_f = k
        barcode_r = reverse_complement(v[0])

        try:
            F_list = l_barcode[barcode_f]
            R_list = r_barcode[barcode_r]
            glue_dic[barcode_f] = [F_list, R_list]

        except KeyError:
            pass
            # print "No pair found for", barcode_f, barcode_r

    return glue_dic


def dic_merge(glue_dic):
    merge_dic={}
    for k in glue_dic.keys():
        v=glue_dic[k]
        merge_dic[k]=[]
        merge_dic[k].append((chain_anchor(v[0])))
        merge_dic[k].append((chain_anchor(v[1])))
    return merge_dic


def merge_filter(merge_dic, covcutoff=3):
    new_dic={}
    for k, v in merge_dic.iteritems():
        if len(v[0])!=1 or len(v[1])!=1:
            v_keep1=[]
            v_keep2=[]
            for bed in v[0]:
                if bed[-1]>=covcutoff:
                    v_keep1.append(bed)
                else:
                    pass
            for bed in v[1]:
                if bed[-1]>=covcutoff:
                    v_keep2.append(bed)
                else:
                    pass
            new_dic[k]=[]
            new_dic[k].append(v_keep1)
            new_dic[k].append(v_keep2)
        else:
            new_dic[k]=v
    return new_dic


def write_tab(merge_dic, out="bac.tab"):
    fw=open(out, "w")
    for k, v in merge_dic.iteritems():
        if len(v[0])==1 and len(v[1])==1:
            chro, start, end ,direction, cov=v[0][0]
            chro_a, start_a, end_a ,direction_a, cov_a=v[1][0]
            #print v[0], v[1]
            v1=[chro, start, end] if direction=="F" else [chro, end, start]
            v2=[chro_a, start_a, end_a] if direction_a=="F" else [chro_a, end_a, start_a]
            v12=v1+v2
            vw=[str(x) for x in v12]
            fw.write("\t".join(vw))
            fw.write("\n")



def flow_bac(barcode_file, bam_f, bam_r, out):
    """
    todo: increase the readability?
    :param barcode_file: "barcodePostpoolce9_BAC_tag.txt"
    :param bam_f: "sF_s.bam"
    :param bam_r: "sR_s.bam"
    :param out: "bac.tab"
    :return:
    """

    tagdict = get_tagpair(barcode_file)

    pos_l=read_position_store(bam_f)
    pos_r=read_position_store(bam_r)

    pos_l_m=read_position_merge(pos_l)
    pos_r_m=read_position_merge(pos_r)

    l_barcode=read_position_barcode(pos_l_m)
    r_barcode=read_position_barcode(pos_r_m)

    glue_dic = glue_FR(tagdict, l_barcode, r_barcode)

    merge_dic=dic_merge(glue_dic)

    new_dic = merge_filter(merge_dic, covcutoff=3)

    write_tab(new_dic, out=out)

    return out



if __name__=="__main__":
    import os
    wkdir="/home/zhaolab1/data/matepair/bac"
    os.chdir(wkdir)
