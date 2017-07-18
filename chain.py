#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/5 16:42
# @Author  : Runsheng     
# @File    : chain.py

"""
The flow for chain contig using longer scaf or chro
(same species assembly or even different species)
"""

import pandas
from operator import itemgetter  # used for sort
from utils import myexe

def read_last_txt(txtfile):
    """
    :param txtfile: the converted lastal file from the last maf alignmnet file
    :return: a dataframe for the alignment
    """
    cmd_rmspace="sed -i \"s/^ //g\" {txtfile}".format(txtfile=txtfile)
    print(cmd_rmspace)
    print(myexe(cmd_rmspace))

    headrow=["score1","chr1", "start1", "len1","direct1", "fulllen1",
             "chr2", "start2", "len2","direct2", "fulllen2","mapping"]
    df = pandas.read_csv(txtfile, sep="\t", comment="#", header=None)
    print(df.shape)
    df.columns = headrow
    print(df.shape)
    return df



def chain_tab(df, chain_interval=15000, anchor_min=200):
    """
    # the most important thing is the output and the input have same format
    # can use both list or both df, but a df combined with a list will be painful
    # judge: add a new chain, or modify chain value
    # rewrite the last line of chain file only to add new boundary
    """
    df = df[(df["len1"] >= anchor_min) & (df["len2"] >= anchor_min)]
    #df = df.sort_values(by="start2")
    # print(df.head)

    chain = df[0:1]

    for n in range(1, len(df)):
        # print list(chain["chr1"])[-1],list(df["chr1"])[n]

        if (list(chain["chr1"])[-1] == list(df["chr1"])[n] and
                    list(chain["chr2"])[-1] == list(df["chr2"])[n] and
                    list(chain["direct1"])[-1] == list(df["direct1"])[n] and
                    list(chain["direct2"])[-1] == list(df["direct2"])[n] and
                        0 <= int(list(df["start1"])[n]) - int(
                            list(chain["start1"])[-1] + list(chain["len1"])[-1]) <= chain_interval and
                        0 <= int(list(df["start2"])[n]) - int(
                            list(chain["start2"])[-1] + list(chain["len2"])[-1]) <= chain_interval):

            chain["len1"][-1:] = list(df["start1"])[n] + list(df["len1"])[n] - list(chain["start1"])[-1]
        else:
            chain = chain.append(df[n:(n + 1)])
    print("The origin alignment have:",df.shape)
    print("The chained alignment have:", chain.shape)
    return chain


if __name__=="__main__":
    import os

    #os.chdir("/home/zhaolab1/myapp/lastchain/test")
    #df=read_last_txt("chro.txt")
    #chain=chain_tab(df)
    #chain.to_csv("chro_chain.txt", sep="\t", header=False, index=False)

    os.chdir("/home/zhaolab1/nanopore/nanoju")
    df=read_last_txt("chro.txt")
    chain=chain_tab(df)
    chain.to_csv("chro_chain.txt", sep="\t", header=False, index=False)

