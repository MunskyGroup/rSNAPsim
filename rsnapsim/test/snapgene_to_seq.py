# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:37:13 2024

@author: willi
"""

import struct

import json
import xmltodict


import html2text

HTML_PARSER = html2text.HTML2Text()
HTML_PARSER.ignore_emphasis = True
HTML_PARSER.ignore_links = True
HTML_PARSER.body_width = 0
HTML_PARSER.single_line_break = True


def parse(val):
    '''parse html'''
    if isinstance(val, str):
        return (HTML_PARSER.handle(val)
                .strip()
                .replace("\n", " ")
                .replace('"', "'"))
    else:
        return val

def get_sequence_from_dnafile(filepath):
    # THIS IS A MODIFIED VERSION OF SNAPGENE READER 
    # https://github.com/IsaacLuo/SnapGeneFileReader/tree/master
    # this only pulls the sequence and name from a given .dna, rSNAPsim assumes 
    # the user knows what they are passing it is an mRNA or CDS or translatable sequence.
    
    f = open(filepath, 'rb')
    
    # read the header first and make sure its snapgene
    
    unpack = lambda size,mode: struct.unpack('>' + mode, f.read(size))[0]
    fb = f.read(1)
    
    if fb != b'\t':
        raise ValueError("Input file is not in SnapGene .dna format")

    
    spacer = unpack(4, 'I')
    title = f.read(8).decode('ascii')
    
    if spacer != 14 or title != 'SnapGene':
        raise ValueError("Input file is not in SnapGene .dna format")

    # features of the snapgene file

    data = dict(is_dna = unpack(2, 'H'),
    exportVersion = unpack(2, 'H'),
    importVersion = unpack(2, 'H'),
    features=[])
    
    bs = []
    while True:
        nb = f.read(1)
        bs.append(nb)
        if nb == b'':
            break
        
        block_size = unpack(4, 'I')
        
        if ord(nb) == 0:
            # read the sequence we still need to pull out the sequence types
            props = unpack(1, 'b')
            data["dna"] = dict(
                topology="circular" if props & 0x01 else "linear",
                strandedness="double" if props & 0x02 > 0 else "single",
                damMethylated=props & 0x04 > 0,
                dcmMethylated=props & 0x08 > 0,
                ecoKIMethylated=props & 0x10 > 0,
                length=block_size - 1
            )
            s = f.read(block_size - 1)
            data["seq"] = s.decode('ascii')     
        else:
            f.read(block_size)
            pass
        
    return data['seq']


filepath = 'C:/Users/willi/Documents/GitHub/rSNAPsim/rsnapsim/test2024/test_gene_files/12xFLAG_ActB_MS2.dna'
bs, data = get_sequence_from_dnafile(filepath)