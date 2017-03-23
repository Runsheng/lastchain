#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/22 14:53
# @Author  : Runsheng     
# @File    : last.py

"""
The wrapper for lastal mapping
"""
from utils import which, myexe
import os

def bin_check(programs):
    """
    check the lastal bins
    :return:
    """
    def bin_is_correct(program):
        if which(program) is None:
            print("No {program} bin in $PATH".format(program=program))
            return False
        else:
            print("{program} bin in $PATH".format(program=program))
            return True

    for program in programs:
        bin_is_correct(program)


def flow_maf(c_target, c_query, wkdir=None, prefix="scafone", core=32, lastal_paral=["-e40", "-150"]):
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
    cmd_lastal="lastal -P{core} {lastal_para} ctarget {c_query} >{prefix}.maf".format(
        core=core, lastal_para=" ".join(lastal_paral), c_query="csp93ctg.fasta", prefix=prefix)
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
    return os.path.join(wkdir, prefix+"text")

def flow(c_target, c_query, wkdir=None):
    """
    contig_reliable: the contig short but accurate, cquery
    contig_long: c_target
    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()

    os.chdir("/home/zhaolab1/myapp/lastchain/test")
    #---
    cmd_lastdb="lastdb ctarget {c_target}".format(c_target="JU.contigs.fasta")
    myexe(cmd_lastdb)

    #---
    cmd_lastal="lastal -P32 -e40 -l50 ctarget {c_query} >scafone.maf".format(c_query="csp93ctg.fasta")
    print myexe(cmd_lastal)

class LastTest(object):
    def __init__(self,strings):
        self.strings=strings

    @staticmethod
    def test_bin_check():
        bin_check(["lastal","lastdb","last-split","maf-sort", "maf-cull","maf-convert"])


    @staticmethod
    def test_flow():
        flow()

if __name__== "__main__":
    pass
    #LastTest.test_bin_check()
    LastTest.test_flow()
