#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/22 14:53
# @Author  : Runsheng     
# @File    : last.py

"""
The wrapper for lastal mapping
"""
from utils import which, myexe, bin_check
import os

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
    cmd_lastdb="lastdb ctarget {c_target}".format(c_target=c_target)
    print(cmd_lastdb)
    myexe(cmd_lastdb)

    #--- use the ctarget as the ref name, evalue
    cmd_lastal="lastal -P{core} {lastal_para} ctarget {c_query} > {prefix}.maf".format(
        core=core, lastal_para=" ".join(lastal_paral), c_query=c_query, prefix=prefix)
    print(cmd_lastal)
    myexe(cmd_lastal)

    # all of these files are temp files
    #---
    cmd_text="""
    last-split {prefix}.maf> {prefix}_s.maf
    maf-sort {prefix}_s.maf > {prefix}_ss.maf
    maf-cull --limit=1 {prefix}_ss.maf >{prefix}_sf.maf
    maf-convert tab {prefix}_sf.maf>{prefix}.txt
    """.format(prefix=prefix)
    print(cmd_text)
    myexe(cmd_text)
    return os.path.join(wkdir, prefix+"_sf.maf")


def reverse_chain(maf, c_target, c_query, kent_path, prefix="scafone", out="reverse", min_score=1000):

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
    wkdir="/home/zhaolab1/nanopore/nanoju"
    flow_maf(c_target="ju_miniasm.fa", c_query="csp93ctg.fasta", wkdir=wkdir, prefix="chro",
             core=40, lastal_paral=["-e40", "-l50"])
