#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/3/22 16:41
# @Author  : Runsheng     
# @File    : scaf.py

"""
The main functions for the pseudo-scaffold

Todo: add the chain function and modify the contig_filter function, for the using of large MB level scafs
"""

# generate a name list from the last output table, select the chr X and length 500
# import pickle as pickle
from collections import OrderedDict
from utils import fasta2dic, reverse_complement, dic2dic, dic2fasta
import os

def contig_filter(filename, mum_cutoff):
    """
    # get the longest MUM for each contig, this MUM can be used to define the location of the contig
    # mum_cutoff is the min length for a MUM record to be retained
    # get a dict to store the contig order

    # add another filter, do not use un or _random contig to layout the contigs,
    # even if they have better mapping than chrs
    """
    contig_d = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            try:
                line_l = line.split("\t")
                chro = line_l[1]
                if (chro == "un") or ("random" in chro) or ("NA" in chro) :
                    pass
                else:
                    name = line_l[6]
                    length = int(line_l[8])

                    if length >= mum_cutoff:
                        if name not in contig_d.keys():
                            contig_d[name] = line_l
                        else:
                            if length > int(contig_d[name][8]):
                                contig_d[name] = line_l
                            else:
                                pass
                    else:
                        pass
            except IndexError:
                pass
    # give some summary info
    print "Total ", len(contig_d), " contigs can be used."
    # if the scaf from the sspace package, will have the length information
    #ks= sorted(contig_d.keys(), key=lambda x: (int(x.split("|")[0].replace("scaffold", ""))))
    #for k in ks:
    #    print contig_d[k]

    return contig_d


def contig_order(filename, mum_cutoff):
    """
    input: the whole genome alignmnet last/chain/net table output for MUMs
    output: the contig order for pseudo scaffolding
    """

    # use pre-define function
    contig_d = contig_filter(filename, mum_cutoff)

    contig_use = []
    for contig_name, contig_l in contig_d.iteritems():
        contig_use.append(contig_l)
    #print(len(contig_use)) # debug only

    # sort according to the chro name and the position
    contig_use_s = sorted(contig_use, key=lambda line: (line[1], int(line[2])))

    chro_d = {}
    for line in contig_use_s:
        chro = line[1]
        contig_name = line[6]
        if chro not in chro_d.keys():
            chro_d[chro] = [contig_name]
        else:
            chro_d[chro].append(contig_name)

    # give some summary info
    for k, v in chro_d.iteritems():
        number = len(v)
        print "There is %d contigs for %s." % (number, k)
    return chro_d, contig_d


def contig_layout(filename, mum_cutoff, ref_dict, prefix, exlist, n100=1):
    """
    lay out the contigs and add 100 "N" between contigs
    generate a gff file to annotate the position of contigs in the pseudo-chromosome
    gff format
    cniX	pseudoscaffold	Contig	2419108	2419128	42	.	.	hid=trf; hstart=1; hend=21
    gff file should be 1 based

    filename="reverse.txt", mum_cutoff=500, ref_dict=fasta2dic("sp9all.fa"), prefix="cni", exlist=[])

    """
    chro_d, contig_d = contig_order(filename=filename, mum_cutoff=mum_cutoff)
    n_100 = "N" * 100 * n100

    chro_gff = []

    chro_allseq = OrderedDict()
    for chro, contigs in chro_d.iteritems():
        chro_seq = []

        start = 1
        end = 1

        for n, contig in enumerate(contigs):  # modified this line to get the last one and hinder the n_100
            is_reverse = contig_d[contig][9]
            if is_reverse == "+":
                contig_seq = ref_dict[contig]
            elif is_reverse == "-":
                contig_seq = reverse_complement(ref_dict[contig])
            chro_seq.append(contig_seq)

            if n != len(contigs) - 1:  # if not last one
                chro_seq.append(n_100)

            end = start + len(contig_seq) - 1  # NOTE: gff is 1 based NOT 0 based
            chrogff_l = [(prefix + chro), "pseudoscaffold", "Contig", str(start), str(end), ".", is_reverse, ".",
                         ("ID=" + contig)]
            chrogff_str = "\t".join(chrogff_l)
            # print(chrogff_str)

            chro_gff.append(chrogff_str)

            start = start + len(contig_seq) + len(n_100) # this is right, equals to start=end +100 +1

        chro_allseq[(prefix + chro)] = "".join(chro_seq)

    # layout un, in unplaced format
    for contig, contig_seq in ref_dict.iteritems():
        if contig not in contig_d.keys() and contig not in exlist:
            chro_allseq[contig] = contig_seq

    # give some summary info
    for chro, seq in chro_allseq.iteritems():
        print chro, len(seq)
    return chro_allseq, chro_gff


def gff2file(chro_gff, out="contig.gff"):
    with open(out,"w") as f:
        for record in chro_gff:
            f.write(record)
            f.write("\n")


def flow_pseudo_scaf(fasta_file, chain_tab, out_prefix="scaf", prefix="", exlist=[], wkdir=None, n100=1):
    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)
    ref_dict=dic2dic(fasta2dic(fasta_file))
    chro_allseq, chro_gff = contig_layout(filename=chain_tab, mum_cutoff=500, ref_dict=ref_dict, prefix=prefix,
                                          exlist=exlist, n100=n100)

    dic2fasta(chro_allseq, "{out_prefix}.fasta".format(out_prefix=out_prefix))
    gff2file(chro_gff, "{out_prefix}.gff".format(out_prefix=out_prefix))


if __name__ == "__main__":

    #wkdir="/home/zhaolab1/myapp/lastchain/test"
    #flow_pseudo_scaf(fasta_file="cni_nanobac.fasta", chain_tab="chro.txt", out_prefix="chro", wkdir=wkdir)

    wkdir="/home/zhaolab1/nanopore/nanoju"
    flow_pseudo_scaf(fasta_file="csp93ctg.fasta", chain_tab="chro_chain.txt", out_prefix="chro", wkdir=wkdir)
