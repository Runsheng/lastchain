#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/7/24 12:14
# @Author  : Runsheng     
# @File    : longscaf.py


"""
Using anchors to reorder the long scaf
"""

from bac import chain_anchor

def chain_anchor(anchors, chain_interval=1000, anchor_min=50):
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