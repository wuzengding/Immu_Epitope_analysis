##########################################################################
# File Name:        4.1.Summary.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Wed 18 May 2022 03:16:02 PM CST
##########################################################################

from argparse import ArgumentParser
import time
import pandas as pd


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-a","--antigenpresent",dest="antigenpresent",
        help="antigen presentation file", required=True)
    parser.add_argument("-i","--immuopeptide",dest="immuopeptide",
        help="immunogeneity peptide file",required=True)
    parser.add_argument("-o","--outfile",dest="outfile",
        help="output file that written result", required=True)
    
    args = parser.parse_args()
    antigenpresent = args.antigenpresent
    immupeptide =  args.immuopeptide
    outfile = args.outfile
    
    antigenpresentdf = pd.read_csv(antigenpresent,names=["seqid","pos","peptide","HLA_subtype","HLA","presentation"])
    immupeptidedf = pd.read_csv(immupeptide,names=["seqid","pos","peptide","HLA_subtype","HLA","immunogeneity"])
    
    mergedf = antigenpresentdf.merge(immupeptidedf,on=["seqid","pos","peptide","HLA_subtype",])
    mergedf.to_csv(outfile)
        
    