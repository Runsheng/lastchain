#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/22 14:53
# @Author  : Runsheng     
# @File    : last.py

"""
The wrapper for lastal mapping
"""
from Utils import which, myexe, bin_check
import os

def flow_maf(c_target, c_query, wkdir=None, prefix="scafone", core=32, lastal_paral=["-e40", "-l50"]):
    """
    try to get a 1:1 alignment, reverse the maf file to filter out all 1:N and N:1 alignemnt

    contig_reliable: the contig short but accurate, cquery
    contig_long: c_target
    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()

    os.chdir(wkdir)
    #---, ctarget is internal reference name
    cmd_lastdb="lastdb ctarget {c_target}".format(c_target=c_target)
    print(cmd_lastdb)
    myexe(cmd_lastdb)

    #--- use the ctarget as the ref name, evalue
    cmd_lastal="lastal -P {core} {lastal_para} ctarget {c_query} > {prefix}.maf".format(
        core=core, lastal_para=" ".join(lastal_paral), c_query=c_query, prefix=prefix)
    print(cmd_lastal)
    myexe(cmd_lastal)

    # all of these files are temp files
    #---
    cmd_text="""
    last-split -m1 {prefix}.maf> {prefix}_s.maf
    maf-sort {prefix}_s.maf > {prefix}_ss.maf
    maf-cull --limit=1 {prefix}_ss.maf >{prefix}_sf.maf
    maf-swap {prefix}_sf.maf | last-split -m1 > reverse.maf
    maf-convert tab reverse.maf>{prefix}.txt
    """.format(prefix=prefix)
    print(cmd_text)
    myexe(cmd_text)
    return os.path.join(wkdir, prefix+".txt")


def __reverse_chain(maf, c_target, c_query, kent_path, prefix="scafone",
                    out="reverse", min_score=1000):

    cmd="""
    export PATH={kent_path}:$PATH
    #convert the maf file to psl file
    #last-split {prefix}.maf>{prefix}_s.maf
    #maf-sort {prefix}_s.maf  >{prefix}_ss.maf
    #maf-cull --limit=1 {prefix}_ss.maf >{prefix}_sf.maf

    maf-convert psl  {maf} > {prefix}.psl

    #prepare the reference to 2-bit format, in two files,  so one can use axtChain correctly
    #in ./raw dict
    faSize {c_target} -detailed > target.sizes
    faSize {c_query} -detailed > query.sizes
    faToTwoBit {c_target} target.fa.2bit
    faToTwoBit {c_query} query.fa.2bit

    #chain order:

    axtChain -psl {prefix}.psl target.fa.2bit query.fa.2bit {prefix}.chain -minScore={min_score} -linearGap=loose -verbose=0
    chainMergeSort {prefix}.chain > {prefix}_s.chain
    chainPreNet {prefix}_s.chain target.sizes query.sizes {prefix}_pre.chain

    #net
    chainNet  {prefix}_pre.chain target.sizes query.sizes target.net query.net
    chainNet  {prefix}_pre.chain target.sizes query.sizes stdout /dev/null |netSyntenic stdin noClass.net

    #nettoaxt
    netToAxt noClass.net {prefix}_pre.chain target.fa.2bit query.fa.2bit stdout | axtSort stdin {out}.axt
    axtToMaf {out}.axt target.sizes query.sizes {out}.maf

    maf-convert tab {out}.maf>{out}.txt
    """.format(prefix=prefix, min_score=min_score, maf=maf,
               kent_path=kent_path, out=out, c_target=c_target, c_query=c_query)

    print(cmd)
    myexe(cmd)
    return out+".txt"


def read_last_to_bed12(lastfile):
    """
    read last txt alignment file, change the file into a two bed combined file

    chr1, start1, end1, chr2, start2, end2, direct2

    chr1 is always forward

    :param lastfile:
    :return:
    """
    f = open(lastfile, "r")
    bed12 = []

    for line in f.readlines():
        if line[0] == "#":
            pass
        else:
            line_l = line.strip().split("\t")
            score1, chr1, start1, len1, direct1, fulllen1, chr2, start2, len2, direct2, fulllen2, mapping, mismap = line_l

            end1 = str(int(start1) + int(len1))
            end2 = str(int(start2) + int(len2))

            if direct2 == "-":
                start2 = str(int(fulllen2) - int(start2) - int(len2))
                end2 = str(int(start2) + int(len2))

            bed_line = (chr1, start1, end1, chr2, start2, end2, direct2)
            bed12.append(bed_line)
    return bed12


def bed2txt(bed12, out="11.txt"):
    with open(out, "w") as fw:
        for line in bed12:
            fw.write("\t".join(line))
            fw.write("\n")



class LastTest(object):
    def __init__(self,strings):
        self.strings=strings

    @staticmethod
    def test_bin_check():
        bin_check(["lastal","lastdb","last-split","maf-sort", "maf-cull","maf-convert"])


    @staticmethod
    def test_flow():
        flow_maf()




if __name__== "__main__":
    #LastTest.test_bin_check()
    #LastTest.test_flow()
    def test1():
        wkdir="/home/zhaolab1/data/matepair/out_hmb/cbrtwo/rc"
        os.chdir(wkdir)
        flow_maf(c_target="../ref/cb4.fasta", c_query="../ref/cni_rc.fa", wkdir=wkdir, prefix="sp9_rc",
              core=30, lastal_paral=["-e40", "-l50"])
        bed12=read_last_to_bed12("sp9_rc.txt")
        bed2txt(bed12, "sp9_rc.bed")

    def reverse2():
        wkdir="/home/zhaolab1/data/matepair/out_hmb/cbrtwo/ref/liftover"
        flow_maf(c_target="csp94.fa", c_query="cni_chro.fa", wkdir=wkdir, prefix="sp9",
                 core=30, lastal_paral=["-e40", "-l50"])
        os.chdir(wkdir)

        kent_path="/home/zhaolab1/app/jk/kentUtils/bin/"
        __reverse_chain(maf="reverse.maf", c_target="csp94.fa", c_query="cni_chro.fa",
                        kent_path=kent_path, out="re", min_score=3000)


