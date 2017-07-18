#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/22 18:58
# @Author  : Runsheng     
# @File    : liftover.py


"""
Using UCSC lift over protocol to transfer the gff file from the origin to the new
require the chain file and the CrossMap package

T rendering based
"""
from utils import which, myexe
from last import bin_check
import os

def chain_wrapper(maf_file, ref_after, ref_before, prefix=None, wkdir=None, minScore=3000):
    """
    :param maf_file:
    :param ref_after:
    :param ref_before:

    :return: the path to the chain file
    """
    if prefix is None:
        prefix=maf_file.split(".")[0]
    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)


    cmd_chain="""
    #PATH setting
    export PATH="/home/zhaolab1/app/jk/kentUtils/bin/":$PATH

    #Convert the maf file to psl file
    maf-convert psl {maf_file} > {prefix}.psl

    #Prepare the reference to 2-bit format, in two files, so can use axtChain correctly

    faSize {ref_after} -detailed > {ref_after}.sizes
    faSize {ref_before} -detailed > {ref_before}.sizes
    faToTwoBit {ref_after} {ref_after}.2bit
    faToTwoBit {ref_before} {ref_before}.2bit

    #chain order
    axtChain -psl {prefix}.psl {ref_before}.2bit {ref_after}.2bit {prefix}.chain \
        -minScore={minScore} -linearGap=loose -verbose=0
    """.format(maf_file=maf_file, ref_after=ref_after, ref_before=ref_before,
               prefix=prefix,minScore=minScore)


    print(cmd_chain)
    myexe(cmd_chain)

    return os.path.join(wkdir, prefix+"chain")


class LastTest(object):
    def __init__(self,strings):
        self.strings=strings

    @staticmethod
    def test_bin_check():
        bin_check(["faSize","faToTwoBit","axtChain","chainMergeSort","maf-convert"])


