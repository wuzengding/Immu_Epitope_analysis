import time,os
from argparse import ArgumentParser
import pandas as pd
import glob
import subprocess
 
#modified at Oct.17 2022 by wuzengding
#modified for adding a argument used to indicate whether make bed file according parsed results or not


if __name__ == '__main__':
    #data_start = time.time()
    #print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Getting Enzyme Site according enzyme result")
    parser = ArgumentParser()
    parser.add_argument("-i","--res_mhci",dest="res_mhcI",
        help="file of MHCi result")
    parser.add_argument("-u","--res_mhcii",dest="res_mhcII",
        help="file of MHCii result")
    parser.add_argument("-n","--res_inhouse",dest="res_inhouse",
        help="file of inhouse result")
    parser.add_argument("-c","--compare",dest="compare",
        help="True|False compare of inhouse with netMHC or not"),
    parser.add_argument("-b","--bedmake",dest="bedmake",
        help="Y/N if it's Y, the program will make bed file",default="N"),
    parser.add_argument("-o","--outpath", dest="outpath")

    
    args = parser.parse_args()
    print(args)
    res_mhcI = args.res_mhcI
    res_mhcII = args.res_mhcII
    res_inhouse = args.res_inhouse
    outpath = args.outpath
    compare = args.compare
    bedmake = args.bedmake

    def prepend_line(file_name, line):
        """ Insert given string as a new line at the beginning of a file """
        # define name of temporary dummy file
        dummy_file = file_name + '.bak'
        # open original file in read mode and dummy file in write mode
        with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
            # Write given line to the dummy file
            write_obj.write(line + '\n')
            # Read lines from original file one by one and append them to the dummy file
            for line in read_obj:
                write_obj.write(line)
        # remove original file
        os.remove(file_name)
        # Rename dummy file as the original file
        os.rename(dummy_file, file_name)



    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    protein_mhctype_list = []
    if res_mhcI:
        os.makedirs(outpath+"/parsed_res_MHCI")
        head_flag = 0
        with open(res_mhcI,"r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                else:
                    if all([element in line for element in ["Pos","MHC","Peptide","Core"]]):
                        header = ",".join(line.split())
                        head_flag += 1
                    else:
                        if line.startswith("------------------------") and head_flag == 1:
                            continue
                        elif line.startswith("------------------------") and head_flag ==2:
                            head_flag =0
                            res_f.close()
                            continue
                        else:
                            if head_flag == 1:
                                mhc_type,gene_name = line.split()[1],line.split()[10]
                                protein_mhctype_list.append((gene_name,"MHCI",mhc_type))
                                res_f = open(os.path.join(outpath,"parsed_res_MHCI",gene_name+"_MHCI_"+mhc_type+".netMHC.csv"),"w")
                                res_f.write(header+"\n")
                                res_f.write(",".join(line.replace("<="," ").split())+"\n")
                                head_flag += 1
                            elif head_flag == 2:
                                res_f.write(",".join(line.replace("<="," ").split())+"\n")
                                
    if res_mhcII:
        os.makedirs(outpath+"/parsed_res_MHCII")
        head_flag = 0
        with open(res_mhcII,"r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                else:
                    if all([element in line for element in ["Pos","MHC","Peptide","Core"]]):
                        header = ",".join(line.split())
                        head_flag += 1
                    else:
                        if line.startswith("------------------------") and head_flag == 1:
                            continue
                        elif line.startswith("------------------------") and head_flag == 2:
                            head_flag = 0
                            res_f.close()
                            continue
                        else:
                            if head_flag == 1:
                                mhc_type,gene_name = line.split()[1],line.split()[6].replace("|","_")
                                protein_mhctype_list.append((gene_name,"MHCII",mhc_type))
                                res_f = open(os.path.join(outpath,"parsed_res_MHCII",gene_name+"_MHCII_"+mhc_type+".netMHC.csv"),"w")
                                res_f.write(header+"\n")
                                res_f.write(",".join(line.replace("<="," ").split())+"\n")
                                head_flag += 1
                            elif head_flag == 2:
                                res_f.write(",".join(line.replace("<="," ").split())+"\n")
    if res_inhouse:
        os.makedirs(outpath+"/parsed_res_InHouse")
        cmd = """
                cut -f1 -d "_" {res_in}|sort |uniq|while read gene_name
                    do
                        echo ${{gene_name}}
                        grep ${{gene_name}} {res_in}|cut -f3 -d","|sort|uniq|while read mhc_type
                            do 
                                echo ${{mhc_type}}
                                grep ${{gene_name}} {res_in}|grep ${{mhc_type}} > \
                                {outpath}/$(echo ${{gene_name}}|sed "s:|:_:g")_${{mhc_type}}.inhouse.csv
                            done
                    done
              """.format(res_in=res_inhouse,outpath=outpath+"/parsed_res_InHouse")
        print(cmd)
        subprocess.call('/usr/bin/bash -c "$cmd"',shell=True,env={'cmd':cmd})
        
    if res_inhouse and compare=="True":
        res_comp_f = open(outpath+"/compare_result.csv","w")
        res_comp_f.write(",".join(["gene_mhctype","peptide_SB_in_netmhc","peptide_WB_in_netmhc",
                                    "SB_and_WB","peptide_POS_inhouse","overlap_SB","overlap_SB_WB"])+"\n")
        for constitute in protein_mhctype_list:
            print(constitute)
            gene_name,mhc,mhc_type= constitute[0],constitute[1],constitute[2]
            mhc_type = "*".join(list(mhc_type.replace("-","")))
            #print(mhc_type)
            netmhc_res_f = glob.glob(os.path.join(outpath,gene_name+"*"+mhc_type+"*"+"netMHC.csv"))
            inhous_res_f = glob.glob(os.path.join(outpath,gene_name+"*"+mhc_type+"*"+"inhouse.csv"))
            #print(netmhc_res_f)
            #print(inhous_res_f)
            netmhc_res = pd.read_csv(netmhc_res_f[0],sep=",",header=0)
            inhous_res = pd.read_csv(inhous_res_f[0],sep=",",names=["gene","Peptide","mhc_type","val"])
            
            netmhc_res_sb = netmhc_res[netmhc_res["BindLevel"]=="SB"]
            netmhc_res_wb = netmhc_res[netmhc_res["BindLevel"]=="WB"]
            netmhc_res_swb = netmhc_res[(netmhc_res["BindLevel"]=="WB") | (netmhc_res["BindLevel"]=="WB")]
            inhous_res_pos = inhous_res[inhous_res["val"]>0.8]
            intersec_res_sb = pd.merge(inhous_res_pos,netmhc_res_sb,on=["Peptide"],how="inner")
            intersec_res_swb = pd.merge(inhous_res_pos,netmhc_res_swb,on=["Peptide"],how="inner")
            
            netmhc_sb = len(netmhc_res_sb)
            netmhc_wb = len(netmhc_res_wb)
            inhous_pos = len(inhous_res_pos)
            intersec_sb = len(intersec_res_sb)
            intersec_swb = len(intersec_res_swb)
            res_comp_f.write(",".join([gene_name+"_"+mhc+"_"+mhc_type,str(netmhc_sb),str(netmhc_wb),str(netmhc_sb+netmhc_wb),str(inhous_pos),str(intersec_sb),str(intersec_swb)])+"\n")
            #print(netmhc_res.head(10))
            #print(inhous_res.head(10))
            #input("Please press the Enter key to proceed")
        res_comp_f.close()
            
    if bedmake == "Y":
        if res_mhcI:
            all_res_files = glob.glob(outpath+"/parsed_res_MHCI/*.netMHC.csv")
            df = pd.concat((pd.read_csv(f) for f in all_res_files),ignore_index=True)
            df.to_csv(outpath+"/parsed_res_MHCI/merge.netMHCI.csv",index=False)
            dfSW = df[(df["BindLevel"]=="SB") | (df["BindLevel"]=="WB")]
            dfSW.loc[:,"protein"] = dfSW["Identity"]
            dfSW.loc[:,"protStart"] = dfSW["Pos"]-1
            dfSW.loc[:,"pepLen"] = dfSW["Peptide"].apply(lambda x:len(x))
            dfSW.loc[:,"protEnd"] = dfSW["Pos"]+dfSW["pepLen"]-1
            dfSW.loc[:,"name"] = dfSW["MHC"]+"_"+df["BindLevel"]
            dfSW.loc[:,"score"] = dfSW["Score_EL"]*1000
            dfSW.loc[:,"strand"] = "."
            dfSW.loc[:,"thickStart"] = dfSW["protStart"]
            dfSW.loc[:,"thickEnd"] = dfSW["protEnd"]
            dfSW.loc[:,"itemRgb"] = "255,0,0"
            dfSW["score"] = dfSW["score"].apply(lambda x: int(x))
            dfSW[["protein","protStart","protEnd","name","score","strand","thickStart","thickEnd","itemRgb"]].to_csv(outpath+"/parsed_res_MHCI/Immunogenicity.netMHCI.bed",sep="\t",index=False,header=False)
            bedheader = '''track name="Immunogenicity_MHCI" description="Immunogenicity prenseted by MHCI" visibility=2 itemRgb="On" useScore=1'''
            prepend_line(outpath+"/parsed_res_MHCI/Immunogenicity.netMHCI.bed",bedheader)
        
        if res_mhcII:
            all_res_files = glob.glob(outpath+"/parsed_res_MHCII/*.netMHC.csv")
            df = pd.concat((pd.read_csv(f) for f in all_res_files),ignore_index=True)
            df.to_csv(outpath+"/parsed_res_MHCII/merge.netMHCII.csv",index=False)
            dfSW = df[(df["BindLevel"]=="SB") | (df["BindLevel"]=="WB")]
            dfSW.loc[:,"protein"] = dfSW["Identity"]
            dfSW.loc[:,"protStart"] = dfSW["Pos"]-1
            dfSW.loc[:,"pepLen"] = dfSW["Peptide"].apply(lambda x:len(x))
            dfSW.loc[:,"protEnd"] = dfSW["Pos"]+dfSW["pepLen"]-1
            dfSW.loc[:,"name"] = dfSW["MHC"]+"_"+df["BindLevel"]
            dfSW.loc[:,"score"] = dfSW["Score_EL"]*1000
            dfSW.loc[:,"strand"] = "."
            dfSW.loc[:,"thickStart"] = dfSW["protStart"]
            dfSW.loc[:,"thickEnd"] = dfSW["protEnd"]
            dfSW.loc[:,"itemRgb"] = "255,0,0"
            dfSW["score"] = dfSW["score"].apply(lambda x: int(x))
            dfSW[["protein","protStart","protEnd","name","score","strand","thickStart","thickEnd","itemRgb"]].to_csv(outpath+"/parsed_res_MHCII/Immunogenicity.netMHCII.bed",sep="\t",index=False,header=False)
            bedheader = '''track name="Immunogenicity_MHCI" description="Immunogenicity prenseted by MHCI" visibility=2 itemRgb="On" useScore=1'''
            prepend_line(outpath+"/parsed_res_MHCII/Immunogenicity.netMHCII.bed",bedheader)