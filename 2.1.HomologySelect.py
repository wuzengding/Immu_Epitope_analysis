##########################################################################
# File Name:        HomologySelect.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Thu 12 May 2022 11:12:25 AM CST
##########################################################################


import time
import os, os.path
from argparse import ArgumentParser
import subprocess

if __name__ == '__main__':
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Homology Selecting .....")
    parser = ArgumentParser()
    parser.add_argument("-f","--blastf",dest="blastf",
        help="file of bastp result",required=True)
    parser.add_argument("-q","--queryseq",dest="queryseq",
        help="digested peptide seq id",required=True)
    parser.add_argument("-m","--mhchlatype",dest="mhchlatype",
        help="HLA type,example:H2-IAb@MHCI",required=True)
    parser.add_argument("-i","--identity",dest="identity",
        help="percent of identity",required=True)
    parser.add_argument("-e","--evalue",dest="evalue",
        help="evalue cutoff",required=False)
    parser.add_argument("-b","-bitscore",dest="bitscore",
        help="bit score cutoff",required=False)
    parser.add_argument("-o","--outfile",dest="outfile",
        help="output file",required=True)
    parser.add_argument("-d","--outdir",dest="outdir",
        help="out path",required=True)
    parser.add_argument("-b","--bedmake",dest="bedmake",
        help="Y/N if it's Y, the program will make bed file",default="N"),
    parser.add_argument("-t","--blastdb",dest="blastdb",
        help="blast database",required=True)

    args = parser.parse_args()
    blastf = args.blastf
    queryseq = args.queryseq
    identity = args.identity
    evalue =  args.evalue
    bitscore = args.bitscore
    outfile = args.outfile
    outdir = args.outdir
    mhchlatype = args.mhchlatype
    blastdb = args.blastdb
    
    #HLA_type,HLA_subtype = mhchlatype.split("@")[0],mhchlatype.split("@")[1]
    #mouseHLAdict={"MHCI":["H2-IAb","H2-IAd"],"MHCII":["H2-Kb","H2-Db","H2-Kd"]}
    
    
    ## get all homoglous seq ids in bast result 
    cmd1 = ''' 
            cat /dev/null > {outdir}/homologous_seqid_homo
            cat {file}|cut -f1 -d'%'|sort|uniq|while read seqid
            do
                #echo ${{seqid}}
                cat /dev/null > {outdir}/temp.homoname
                grep "${{seqid}}" {file}|cut -f2 -d'%' |sort|uniq | tr ':' '\n' >> {outdir}/temp.homoname
                #cat /mnt/data2/wuzengding/00.database/index_mus_blastdb/tr_mouse_Uncharacterized_protein.fasta >> {outdir}/temp.homoname
                awk 'NR==FNR{{a[$1];next}} !($2 in a)' {outdir}/temp.homoname {file} | grep "${{seqid}}" | awk -F'\t' '{{if ($3 >{identity}) {{print $1}}}}'|sort |uniq>> {outdir}/homologous_seqid_homo
            done
           '''
    print(cmd1.format(file=blastf,identity=identity,outdir=outdir))
    os.system(cmd1.format(file=blastf,identity=identity,outdir=outdir))
    
    
    ##get all seq ids in query file
    cmd2 = ''' cat /dev/null > {outdir}/homologous_seqid_query
               cat {queryseq}|grep ">" |cut -d'>' -f2 |sort|uniq >> {outdir}/homologous_seqid_query
               cat {outdir}/homologous_seqid_homo {outdir}/homologous_seqid_query > {outdir}/homologous_seqid '''
    os.system(cmd2.format(queryseq=queryseq,outdir=outdir))
    
    os.system('''cat /dev/null > {outfile}'''.format(outfile=outfile))

    cmd3_temp = '''cat {outdir}/homologous_seqid|sort| uniq -u |while read seqid
              do
                seqseq=$(grep -A1 $seqid {queryseq}|grep -v $seqid)
                #echo $seqseq
                if [ ${{#seqseq}} -ge 8 ] && [ ${{#seqseq}} -le 12 ];then
                    #echo "haha"
                    if [ {mhchlatype} == "mouse" ];then
                        #echo "yes"
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-Kb'\t'MHCI >> {outfile}
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-Db'\t'MHCI >> {outfile}
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-Kd'\t'MHCI >> {outfile}
                    else
                        echo ${{seqid}}'\t'${{seqseq}}'\t'${{{mhchlatype}%@*}}'\t'${{{mhchlatype}#*@}} >> {outfile}
                    fi
                fi
                
                if [ ${{#seqseq}} -ge 9 ] && [ ${{#seqseq}} -le 30 ];then
                    if [ {mhchlatype} == "mouse" ];then
                        #echo "yes1"
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-IAb'\t'MHCII >> {outfile}
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-IAd'\t'MHCII >> {outfile}
                    else
                        echo ${{seqid}}'\t'${{seqseq}}'\t'${{{mhchlatype}%@*}}'\t'${{{mhchlatype}#*@}} >> {outfile}
                    fi
                fi
                
              done
                '''.format(outdir=outdir,queryseq=queryseq,outfile=outfile,mhchlatype=mhchlatype)
    cmd3 = '''cat {outdir}/homologous_seqid|sort| uniq -u |while read seqid
              do
                seqseq=$(grep -A1 $seqid {queryseq}|grep -v $seqid)
                #echo $seqseq
                if [ ${{#seqseq}} -ge 8 ] && [ ${{#seqseq}} -le 12 ];then
                    #echo "haha"
                    if [ {mhchlatype} == "mouse" ];then
                        #echo "yes"
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-Kb'\t'MHCI >> {outfile}
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-Db'\t'MHCI >> {outfile}
                    else
                        echo ${{seqid}}'\t'${{seqseq}}'\t'${{{mhchlatype}%@*}}'\t'${{{mhchlatype}#*@}} >> {outfile}
                    fi
                fi
                
                if [ ${{#seqseq}} -ge 9 ] && [ ${{#seqseq}} -le 30 ];then
                    if [ {mhchlatype} == "mouse" ];then
                        #echo "yes1"
                        echo ${{seqid}}'\t'${{seqseq}}'\t'H2-IAb'\t'MHCII >> {outfile}
                    else
                        echo ${{seqid}}'\t'${{seqseq}}'\t'${{{mhchlatype}%@*}}'\t'${{{mhchlatype}#*@}} >> {outfile}
                    fi
                fi
                
              done
                '''.format(outdir=outdir,queryseq=queryseq,outfile=outfile,mhchlatype=mhchlatype)
    #print(cmd3)
    subprocess.call('/usr/bin/bash -c "$cmd"',shell=True,env={'cmd':cmd3})
    