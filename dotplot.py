#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/17 12:31
# @Author  : Runsheng     
# @File    : dotplot.py

"""
combine two chromosome sets into two continuous single contigs, and make alignment and draw figs
"""

from utils import fasta2dic, myexe
from last import flow_maf
import os

def combine_reads(record_dict, tree_keys=None, N_no=2000, prefix="combine"):
    fasta = []
    bed = []

    start = 1
    if tree_keys is None:
        tree_keys = sorted(record_dict.keys())

    for name in tree_keys:
        seq = record_dict[name]
        fasta.append(str(seq.seq))

        if len(bed) == 0:
            end = start + len(seq)
        else:
            start = end + N_no
            end = start + len(seq)
        bed4 = (prefix, str(start), str(end), name)
        bed.append(bed4)
    Ns = "N" * N_no
    fasta_n = Ns.join(fasta)

    return fasta_n, bed


def write_bed(bed, prefix="combine"):
    filename=prefix+".bed"
    with open(filename, "w") as fw:
        for i in bed:
            fw.write("\t".join(i))
            fw.write("\n")
    return filename


def write_combine_fasta(fasta_n, prefix="combine"):
    filename=prefix+".fasta"
    with open(filename, "w") as fw:
        fw.write(">"+prefix+"\n")
        fw.write(fasta_n)
    return filename


def flow_one(record_dict, prefix="combine", tree_keys=None):
    """
    flow to combine and wirte fasta and bed4 files
    :return: filename for fasta and bed
    """
    fasta_n, bed=combine_reads(record_dict=record_dict,tree_keys=tree_keys, prefix=prefix)
    bed_file=write_bed(bed, prefix=prefix)
    fasta_file=write_combine_fasta(fasta_n, prefix=prefix)

    return fasta_file, bed_file


def flow_maf(c_target, c_query, wkdir=None, prefix="scafone", core=32, lastal_paral=["-e40", "-l50"]):
    """
    contig_reliable: the contig short but accurate, cquery
    contig_long: c_target
    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()

    os.chdir(wkdir)
    #---, ctarget is internal reference name
    cmd_lastdb="lastdb -P{core} -uMAM4 -R01 ctarget {c_target}".format(c_target=c_target,core=core)
    print(cmd_lastdb)
    myexe(cmd_lastdb)

    #--- train lastal, new in versin 8
    cmd_lasttrain="last-train -P{core} --revsym --matsym --gapsym -E0.05 ctarget {c_query} > train.mat".format(
        c_query=c_query, core=core)
    print(cmd_lasttrain)
    myexe(cmd_lasttrain)
    #--- use the ctarget as the ref name, evalue
    cmd_lastal="lastal -P{core} {lastal_para} -p train.mat ctarget {c_query} > {prefix}.maf".format(
        core=core, lastal_para=" ".join(lastal_paral), c_query=c_query, prefix=prefix)
    print(cmd_lastal)
    myexe(cmd_lastal)

    # all of these files are temp files
    #---
    cmd_text="""
    maf-sort {prefix}.maf > {prefix}_ss.maf
    last-split -m1 {prefix}_ss.maf >{prefix}_sf.maf
    maf-convert tab {prefix}_sf.maf>{prefix}.txt
    sed -i "s/^ //g" {prefix}.txt
    """.format(prefix=prefix)
    print(cmd_text)
    myexe(cmd_text)
    return os.path.join(wkdir, prefix+".txt")



def flow_dotplot(c_target, c_query, wkdir=None, key_target=None, key_query=None):
    """
    Combine the target and query genome as one, marked in bed file, and draw the dotplot
    :param c_target:
    :param c_query:
    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)

    target_dict=fasta2dic(c_target)
    query_dict=fasta2dic(c_query)

    target_fasta, target_bed=flow_one(target_dict, prefix="target", tree_keys=key_target)
    query_fasta, query_bed=flow_one(query_dict, prefix="query",tree_keys=key_query)

    flow_maf(c_target=target_fasta, c_query=query_fasta, prefix="1vs1",core=40, lastal_paral=["-E0.05", "-m100", "-C2"])

if __name__=="__main__":
    c_target="cb4.fasta"
    c_query="pseudo.fasta"
    wkdir="/home/zhaolab1/data/matepair/pseudo"
    keys=["I", "II", "III", "IV", "V", "X"]

    flow_dotplot(c_target=c_target, c_query=c_query, wkdir=wkdir,key_target=keys, key_query=keys)


