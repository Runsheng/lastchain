#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/17 11:16
# @Author  : Runsheng     
# @File    : flow.py


from last import flow_maf, reverse_chain
from chain import read_last_txt, chain_tab
from scaf import flow_pseudo_scaf

import os

def flow_last_chain_pseudoscaf(c_target, c_query, kent_path, wkdir=None, prefix="pseudo",
                               core=40, lastal_para=["-e40", "-l100"]):
    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)

    maf=flow_maf(c_target, c_query, wkdir, prefix, core, lastal_para)

    maf_txt=reverse_chain(maf=maf, kent_path=kent_path, c_target=c_target, c_query=c_query,
                          prefix=prefix, out="reverse",min_score=5000)

    df=read_last_txt(maf_txt)
    chain=chain_tab(df)
    chain_txt=prefix+"_chain.txt"
    chain.to_csv(chain_txt, sep="\t", header=False, index=False)
    flow_pseudo_scaf(fasta_file=c_query, chain_tab=chain_txt, out_prefix=prefix,wkdir=wkdir)


if __name__=="__main__":
    wkdir="/home/zhaolab1/data/matepair/pseudo"
    kent_path="/home/zhaolab1/app/jk/kentUtils/bin/"
    c_target="cb4.fasta"
    c_query="cni_nano_mb.fa"
    prefix="pseudo"
    os.chdir(wkdir)
    #flow_pseudo_scaf(fasta_file=c_query, chain_tab="pseudo_chain.txt", out_prefix=prefix,wkdir=wkdir, n100=10)
    # flow_last_chain_pseudoscaf(c_target=c_target, c_query=c_query,wkdir=wkdir, kent_path=kent_path)
    maf_txt = reverse_chain(maf="pseudo_sf.maf", kent_path=kent_path, c_target=c_target, c_query=c_query,
                            prefix="first", out="reverse", min_score=1000)

    df = read_last_txt("reverse.txt")
    chain = chain_tab(df)
    chain_txt = prefix + "_chain.txt"
    chain.to_csv(chain_txt, sep="\t", header=False, index=False)
    flow_pseudo_scaf(fasta_file=c_query, chain_tab=chain_txt, out_prefix=prefix, wkdir=wkdir)


    #c_query="pseudo.fasta"
    #flow_maf(c_target,c_query,wkdir, prefix="chro",lastal_paral=["-e40", "-l500"])