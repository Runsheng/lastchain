#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/5 16:42
# @Author  : Runsheng     
# @File    : chain.py

"""
The flow for chain contig using longer scaf or chro
(same species assembly or even different species)

TODO: re-write the chain anchor function
"""

import pandas
from operator import itemgetter  # used for sort
from Utils import myexe
from operator import itemgetter


def read_last_txt(txtfile):
    """
    :param txtfile: the converted lastal file from the last maf alignmnet file
    :return: a dataframe for the alignment
    """
    cmd_rmspace="sed -i \"s/^ //g\" {txtfile}".format(txtfile=txtfile)
    print(cmd_rmspace)
    print(myexe(cmd_rmspace))

    headrow=["score1","chr1", "start1", "len1","direct1", "fulllen1",
             "chr2", "start2", "len2","direct2", "fulllen2","mapping", "mismap"]
    df = pandas.read_csv(txtfile, sep="\t", comment="#", header=None)
    print(df.shape)
    if df.shape[1]==len(headrow):
        df.columns = headrow
    else:
        print("Length unmatch! Please check the txt alignment file!")
    return df


def read_blast_m6(txtfile):
    """
    :param textfile : the blast m6 file
    :return: a dataframe for the alignment
    """


def df2anchor(df):
    """

    :param df:
    :return: anchor dict for each scaf, as {chr2:[(chr2, chr1, start1, end1, direction)...]...}
    """

    """
    get a dict for each scaf, with all other chrom/scaf
    """
    anchor_d={}

    for index, row in df.iterrows():
        chr2=row[6]
        chr1=row[1]
        start1=int(row[2])
        len1=int(row[3])
        end1=start1+len1
        direct1=row[4]
        direct2=row[9]
        direction="F" if direct1==direct2 else "R"
        anchor_line=(chr2, chr1, start1, end1, direction)

        try:
            anchor_d[chr2].append(anchor_line)
        except KeyError:
            anchor_d[chr2]=[]
            anchor_d.append(anchor_line)

    return anchor_d


def chain_anchor(anchors, chain_interval=5000, anchor_min=50):
    """
        Input anchor is the list with 5 element tuples with sorted
        'C00126:229:C8T9TANXX:8:1216:8219:61822BARCODECCCATGGCTGTTGGAGACAG','ctg100000399369',59929,60239,'R'),

        return: 'ctg100000399369',59929,60239,'R', count
    """
    anchor = sorted(anchors, key=itemgetter(1, -1, 2, 3))
    chains = []
    count = 0
    for n, line in enumerate(anchor):
        if n == 0:
            read, chro, start, end, direction = line
            chain_start = start
            chain_end = end
            chain_direction = direction
            # print "first chro", chro
            count += 1
        elif n > 0:
            read, chro, start, end, direction = line
            read_p, chro_p, start_p, end_p, direction_p = anchor[n - 1]  # previous one

            if chro == chro_p and direction == direction_p:
                if end_p >= chain_end and end_p - chain_end <= chain_interval:
                    chain_end = end_p
                    count += 1
                elif end_p > chain_end and end_p - chain_end > chain_interval:
                    chains.append((chro, chain_start, chain_end, chain_direction,
                                   count))  # if break, add the previous chain to list
                    # print "New chr", chro
                    chain_start = start
                    chain_end = end
                    chain_direction = direction
                    count = 1

            elif chro != chro_p or direction != direction_p:
                chains.append(
                    (chro, chain_start, chain_end, chain_direction, count))  # if break, add the previous chain to list
                # print "New chr", chro
                chain_start = start
                chain_end = end
                chain_direction = direction
                count = 1

    chains.append((chro, chain_start, chain_end, chain_direction, count))  # if not break, add the final one into list

    return chains


def __chain_tab(df, chain_interval=15000, anchor_min=200):
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


def chain_long():
    pass


def read_last_txt_im(txtfile):
    pass



if __name__=="__main__":
    import os

    #os.chdir("/home/zhaolab1/myapp/lastchain/test")
    #df=read_last_txt("chro.txt")
    #chain=chain_tab(df)
    #chain.to_csv("chro_chain.txt", sep="\t", header=False, index=False)

    os.chdir("/home/zhaolab1/nanopore/nanoju")
    df=read_last_txt("chro.txt")
    #chain=chain_tab(df)
    #chain.to_csv("chro_chain.txt", sep="\t", header=False, index=False)

