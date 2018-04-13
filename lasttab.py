#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/8/8 12:14
# @Author  : Runsheng     
# @File    : lasttab.py


"""
The class used to parse the last tab file, store and merge all needed field

"""
from collections import namedtuple
from operator import attrgetter

last_fields=["score1","chr1", "start1", "len1","direct1", "fulllen1",
         "chr2", "start2", "len2","direct2", "fulllen2","mapping", "mismap"]
last_record=namedtuple("lasttab", last_fields)


class LastTAB(object):
    """
    last tab file parser
    """

    def __init__(self, tabfile):
        self.name=tabfile
        self._get_atr()


    def _get_atr(self):
        """
        pare: a last tab file
        get a nametuple of all lines in self.atr
        """
        self.atr=[]
        with open(self.name, "r") as f:
            for line in f.readlines():
                if line[0]=="#" or line[0]==" ":
                    pass
                else:
                    line_t=line.strip().split("\t")
                    last_atr=last_record(*line_t)
                    self.atr.append(last_atr)

    def get_link(self):
        """
        using self.atr to get a chro:chro link
        :return:
        """
        for tab in self.atr:
            pass





if __name__=="__main__":

    # test code

    import os
    os.chdir("/home/zhaolab1/syn_jcvi/sample/three/genescaf/align")



