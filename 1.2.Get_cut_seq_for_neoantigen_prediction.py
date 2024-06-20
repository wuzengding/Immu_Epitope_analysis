##########################################################################
# File Name:        Get_cut_seq_for_neoantigen_prediction.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 16 May 2022 10:17:31 AM CST
# modifed time:     Thirsday 08 Sept 2022  by wuzengding  
# 修改原因          修改后的脚本可以直接接受fasta序列，生成可以直接用于inohuse 新抗原呈递预测的格式
##########################################################################


import time,os
from argparse import ArgumentParser
import pandas as pd

if __name__ == '__main__':
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Getting Enzyme Site according enzyme result")
    parser = ArgumentParser()
    parser.add_argument("-f","--inputfasta",dest="inputfasta",
        help="fasta file of input")
    parser.add_argument("-o","--outpath",dest="outpath",
        help="file of output which can be as input of inhouse AP prediction software")
    parser.add_argument("-b","--binsizelist",dest="binsizelist",
        help="bin size list, eg: MHCI:9,10_MHCII:13,15,which means MCHI peptide length in 9,10aa, and MHCII peptide length is 13,15aa ",required=True)
    parser.add_argument("-g","--genetypeofMHC",dest="genetypeofMHC",
        help="gene type list, eg: MHCI:H2-Kb,H2-Db_MHCII:H2-IAb")
    parser.add_argument("-s","--stepsize",dest="stepsize",
        help="step size",required=True)
        
    args = parser.parse_args()
    inputfasta = args.inputfasta
    outpath = args.outpath
    stepsize = int(args.stepsize)
    ##解析所有要切出的MHCI 和MHCII的长度参数
    binsizelist = args.binsizelist
    MHCIBinsize = [int(i)  for i in binsizelist.split("_")[0].replace("MHCI:","").split(",")]
    MHCIIBinsize = [int(i) for i in binsizelist.split("_")[1].replace("MHCII:","").split(",")]
    print("peptide for MHCI length  is {MHCIBinsize}".format(MHCIBinsize=MHCIBinsize))
    print("peptide for MHCII length is {MHCIIBinsize}".format(MHCIIBinsize=MHCIIBinsize))
    ##解析所有要分析的MHCI 和MHCII的基因型参数
    genetypeofMHC = args.genetypeofMHC
    MHCItype = [tp for tp in genetypeofMHC.split("_")[0].replace("MHCI:","").split(",")]
    MHCIItype =[tp for tp in genetypeofMHC.split("_")[1].replace("MHCII:","").split(",")]
    print("MHCI genetype is {MHCItype}".format(MHCItype=MHCItype))
    print("MHCII genetype is {MHCIItype}".format(MHCIItype=MHCIItype))



    def cut_sequence(seq,stepsize,binsize):
        """
        return list of seq window slide
        seq must be string
        """
        return [seq[i:i+binsize] for i in range(0, len(seq)-binsize, stepsize)]
    
    def writeout(seq_id,seq_list,outf,MHC,MHCtype):
        """
        输出对应的序列到结果文件中
        """
        for i,sequence in enumerate(seq_list):
            for genetype in MHCtype:
                outf.write("\t".join([seq_id+"_"+str(i),sequence,genetype,MHC])+"\n")
                
                
                
    if len(MHCIBinsize) > 0: 
        """
        如果MHCI 需要有对应输出的话，才会得到MCHI对应序列，
        只输入MHCI的genetype不会有得到对应序列
        """
        for binsize in MHCIBinsize:
            outf_name = "MHCI_{binsize}.csv".format(binsize=binsize)
            outf = open(os.path.join(outpath,outf_name),"w")
            with open(inputfasta,"r") as f:
                seq_seq = ""
                first_seq_tag = 1
                for line in f.readlines():
                    if line.startswith(">"):
                        if first_seq_tag == 1:
                            first_seq_tag = 0
                        else:
                            seq_list = cut_sequence(seq_seq,stepsize,binsize)
                            writeout(seq_id,seq_list,outf,"MHCI",MHCItype)
                        seq_id = line.strip().replace(">","")
                    else:
                        seq_seq += line.strip()
                
                ##输出最后一行序列，最后一行序列容易漏输出
                seq_list = cut_sequence(seq_seq,stepsize,binsize)
                writeout(seq_id,seq_list,outf,"MHCI",MHCItype)
                
    if len(MHCIIBinsize) > 0: 
        """
        如果MHCII 需要有对应输出的话，才会得到MCHI对应序列，
        只输入MHCI的genetype不会有得到对应序列
        """
        for binsize in MHCIIBinsize:
            outf_name = "MHCII_{binsize}.csv".format(binsize=binsize)
            outf = open(os.path.join(outpath,outf_name),"w")
            with open(inputfasta,"r") as f:
                seq_seq = ""
                first_seq_tag = 1
                for line in f.readlines():
                    if line.startswith(">"):
                        if first_seq_tag == 1:
                            first_seq_tag = 0
                        else:
                            seq_list = cut_sequence(seq_seq,stepsize,binsize)
                            writeout(seq_id,seq_list,outf,"MHCII",MHCIItype)
                        seq_id = line.strip().replace(">","")
                    else:
                        seq_seq += line.strip()
                
                ##输出最后一行序列，最后一行序列容易漏输出
                seq_list = cut_sequence(seq_seq,stepsize,binsize)
                writeout(seq_id,seq_list,outf,"MHCII",MHCIItype)
