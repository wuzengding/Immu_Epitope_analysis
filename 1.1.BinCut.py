##########################################################################
# File Name:        BinCut.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Wed 11 May 2022 10:44:52 AM CST
##########################################################################

import os
import time
from argparse import ArgumentParser

if __name__ == '__main__':
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Bin Cutting .....")
    parser = ArgumentParser()
    parser.add_argument("-f","--fasta",dest="fafile",
        help="fasta file that need been cutted by bin step",required=True)
    parser.add_argument("-b","--binsize",dest="binsize",
        help="bin step size",required=True)
    parser.add_argument("-o","--out1",dest="outfile1",
        help="output file that written cutted sequence as table separate format",required=True)
    parser.add_argument("-u","--out2",dest="outfile2",
        help="output file that written cutted sequence as fa format",required=True)
    
    args = parser.parse_args()
    fafile = args.fafile
    binsize = args.binsize
    outfile1 = args.outfile1
    outfile2 = args.outfile2
    binsize = int(binsize)
    step=40
    outf1 = open(outfile1,"w")
    outf2 = open(outfile2,"w")

    def window_generator(seq, window_lenth=16,step=1):
        """
        return list of seq window slide
        seq is string
        """
        return [seq[i:i+win_len] for i in range(0,len(seq),step)]
    
    first_seq=True
    with open(fafile,"r") as f:
        fileline = f.readlines()
        linetotal = len(fileline)
        for iline in range(linetotal):
            #print(iline,linetotal)
            line = fileline[iline]
            if iline == 0:
                seq_id = line.strip().split(">")[1].replace(" ","_").replace(",","-")
                seq_seq = ""
            elif iline < linetotal-1:
                if line.startswith(">"):
                    length = len(seq_seq)
                    if length > binsize:
                        for i in range(0,length,step):
                            seq_name = seq_id+":"+str(i)+"-"+str(i+binsize)
                            seq_seq_step = seq_seq[i:i+binsize]
                            if (i+binsize)>= length:
                                seq_name = seq_id+":"+str(i)+"-"+str(length)
                                seq_seq_step = seq_seq[i:length]
                                outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                                outf2.write("\t".join([seq_name,seq_seq_step])+"\n")
                                break
                            else:
                                outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                                outf2.write("\t".join([seq_name,seq_seq_step])+"\n")
                    else:
                        outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                        outf2.write("\t".join([seq_id,seq_seq])+"\n")
                    seq_id = line.strip().split(">")[1].replace(" ","_").replace(",","-")
                    seq_seq = ""
                else:
                    seq_seq += line.strip()
            elif iline == linetotal-1:
                seq_seq += line.strip()
                length = len(seq_seq)
                if length > binsize:
                    for i in range(0,length,step):
                        seq_name = seq_id+":"+str(i)+"-"+str(i+binsize)
                        seq_seq_step = seq_seq[i:i+binsize]
                        if (i+binsize)>= length:
                            seq_name = seq_id+":"+str(i)+"-"+str(length)
                            seq_seq_step = seq_seq[i:length]
                            outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                            outf2.write("\t".join([seq_name,seq_seq_step])+"\n")
                            break
                        else:
                            #print("\t".join([seq_name,seq_seq_step])+"\n")
                            outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                            outf2.write("\t".join([seq_name,seq_seq_step])+"\n")
                else:
                    outf1.write("\n".join([">"+seq_name,seq_seq_step])+"\n")
                    outf2.write("\t".join([seq_id,seq_seq])+"\n")
    outf1.close()
    outf2.close()