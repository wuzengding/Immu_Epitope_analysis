##########################################################################
# File Name:        GetEnzymeSite.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 16 May 2022 10:17:31 AM CST
##########################################################################

##Version=V1

import time,os
from argparse import ArgumentParser
import pandas as pd

if __name__ == '__main__':
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Getting Enzyme Site according enzyme result")
    parser = ArgumentParser()
    parser.add_argument("-f","--inputfile",dest="inputfile",
        help="file of input")
    parser.add_argument("-o","--outfile",dest="outfile",
        help="file of output")
    parser.add_argument("-b","--binsize",dest="binsize",
        help="bin step size",required=True)
        
    args = parser.parse_args()
    inputfile = args.inputfile
    outfile = args.outfile
    binsize = int(args.binsize)
    
    ## 8aa <= length <= 12aa
    ## 9aa <= length <= 30aa
    MHCI_min,MHCI_max = 8,12
    MHCII_min,MHCII_max = 9,25
    
    cutoff = 0.5
    outf = open(outfile,"w")
    seqid = ""
    wholeseq = ""
    sitelist = [0]
    with open(inputfile,"r") as f:
        fileline = f.readlines()
        linetotal = len(fileline)
        #print(linetotal)
        for iline in range(linetotal):
            print(iline,linetotal)
            line = fileline[iline]
            seq_id,seq_seq,seq_res = line.strip().rsplit(':',1)[0],line.strip().split(',')[1], float(line.strip().split(',')[2])
            #print(seq_id,seq_seq,seq_res)
            if iline == 0:
                seqid = seq_id
                wholeseq = seq_seq
                if seq_res >= cutoff:
                    print("yes")
                    site = binsize//2
                    sitelist.append(site)
            elif iline < linetotal-1:
                if seq_id == seqid:
                    seqid = seq_id
                    wholeseq += seq_seq[-1]
                    if seq_res >= cutoff:
                        print("yes",wholeseq)
                        site = len(wholeseq) - binsize//2
                        sitelist.append(site)
                elif seq_id != seqid:
                    #sitelist.append(len(wholeseq))
                    #print(sitelist)
                    #for i in range(len(sitelist)-2):
                    #    for j in range(i+1,len(sitelist)-2):
                    #        pos1,pos2 = sitelist[i],sitelist[j]
                    #        seqget = wholeseq[pos1:pos2]
                    #        if MHCI_min <=len(seqget) <= MHCII_max:
                    #            outf.write('\n'.join([">"+seqid+":"+str(pos1)+"_"+str(pos2),seqget])+"\n")
                    #        elif len(seqget) > MHCII_max:
                    #            break
                    for i in range(len(wholeseq)):
                        for j in range(MHCI_min,MHCII_max+1):
                            pos1,pos2 = i,i+j
                            seqget = wholeseq[pos1:pos2]
                            outf.write('\n'.join([">"+seqid+":"+str(pos1)+"_"+str(pos2),seqget])+"\n")

                    wholeseq = seq_seq
                    sitelist = [0]
                    seqid = seq_id
                    #print(seqid)
            elif iline == (linetotal-1):
                #sitelist.append(len(wholeseq))
                #print(sitelist)
                #for i in range(len(sitelist)-2):
                #    for j in range(i+1,len(sitelist)-2):
                #        pos1,pos2 = sitelist[i],sitelist[j]
                #        seqget = wholeseq[pos1:pos2]
                #        if MHCI_min <=len(seqget) <= MHCII_max:
                #            outf.write('\n'.join([">"+seqid+"_"+str(pos1)+"_"+str(pos2),seqget])+"\n")
                #        elif len(seqget) > MHCII_max:
                #            break
                
                for i in range(len(wholeseq)):
                    for j in range(MHCI_min,MHCII_max+1):
                        pos1,pos2 = i,i+j
                        seqget = wholeseq[pos1:pos2]
                        outf.write('\n'.join([">"+seqid+":"+str(pos1)+"_"+str(pos2),seqget])+"\n")
    outf.close()
    
            
              
    #enzyme_res = pd.read_csv(input,header=["seq_id","seq_seq","seq_res"],)
    