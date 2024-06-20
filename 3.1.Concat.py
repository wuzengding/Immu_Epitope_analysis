import sys
import requests
from argparse import ArgumentParser
import time


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-r","--ruler",dest="ruler",
        help="Rule name that ask for eg:  \
               EnzymeDigestion| AntigenPresentation| Immunogenicity ", required=True)
    parser.add_argument("-i","--inputfile",dest="inputfile",
        help="input file that need analysis",required=True)
    parser.add_argument("-o","--outfile",dest="outfile",
        help="output file that written result", required=True)
    
    args = parser.parse_args()
    ruler = args.ruler
    inputfile = args.inputfile
    outfile = args.outfile

    url1_18aa2vs1_mhc12_v1 = "https://51120m7l36.zicp.vip:443/enzyme/"
    url1_12aa2vs3_mhc1_v1 = "https://51120m7l36.zicp.vip:443/enzyme0601/"
    url_20220625 = "https://51120m7l36.zicp.vip:443/enzyme/flatten_80/"
    aa80_cleav_9606_I_win80setp40_matrix_mergeCut="https://51120m7l36.zicp.vip:443/enzyme/matrix_allStar_V2/"
    #aa80_cleav_9606_I_win80setp40_matrix_mergeCut = "http://192.168.2.253:8000/enzyme/matrix_allStar_V2/"
    url2 = "https://51120m7l36.zicp.vip:443/antigen_presentation/"
    url3 = "https://51120m7l36.zicp.vip:443/immunogenicity/"
    
    data_start = time.time()
    with open(inputfile,"rb") as f:
        if ruler=="EnzymeDigestion":
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " EnzymeDigestion analysis .....")
            res = requests.post(url=aa80_cleav_9606_I_win80setp40_matrix_mergeCut, files={"fileb": f}, timeout=None)
    
        elif ruler=="AntigenPresentation":
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " AntigenPresentation analysis .....")
            res = requests.post(url=url2, files={"fileb": f}, timeout=None)
        
        elif ruler=="Immunogenicity":
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Immunogenicity analysis .....")
            res = requests.post(url=url3, files={"fileb": f}, timeout=None)
    
        res = res.text
        print(res)
        outf = open(outfile,"w")
        outf.write(res)
        outf.close()
