#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/24 10:25
# @Author  : Runsheng     
# @File    : flow_test.py


from chain import df2anchor, chain_anchor

import os


def flow_chain_anchor():
    pass






if __name__=="__main__":
    wkdir="/home/zhaolab1/data/matepair/pseudo"
    kent_path="/home/zhaolab1/app/jk/kentUtils/bin/"
    c_target="cb4.fasta"
    c_query="cni_nano_mb.fa"
    prefix="pseudo"
    os.chdir(wkdir)

    #flow_pseudo_scaf(fasta_file=c_query, chain_tab="pseudo_chain.txt", out_prefix=prefix,wkdir=wkdir, n100=10)
    #maf_txt = __reverse_chain(maf="pseudo_sf.maf", kent_path=kent_path, c_target=c_target, c_query=c_query,
    #                        prefix="first", out="reverse", min_score=1000)

    df = read_last_txt("reverse.txt")
    chain = __chain_tab(df)
    chain_txt = prefix + "_chain.txt"
    chain.to_csv(chain_txt, sep="\t", header=False, index=False)
    flow_pseudo_scaf(fasta_file=c_query, chain_tab=chain_txt, out_prefix=prefix, wkdir=wkdir)