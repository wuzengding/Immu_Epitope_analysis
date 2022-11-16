##########################################################################
# File Name:        00.vaccine_sequence_design.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Tue 10 May 2022 03:15:55 PM CST
##########################################################################

configfile: os.path.join("/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/Config.pipline.json")
#print(config)
inhousetoolsconcat = config["tools"]["inhousetoolsconcat"]
BinCut = config["tools"]["BinCut"]
EnzymeDigestion = config["tools"]["EnzymeDigestion"]
GetCutSeq = config["tools"]["GetCutSeq"]
Homology = config["tools"]["Homology"]
AntigenPresentation = config["tools"]["AntigenPresentation"]
Immunogenicity = config["tools"]["Immunogenicity"]
Summary = config["tools"]["summary"]
#Population = config["tools"]["Population"]
RawPeptide = config["RawPeptide"]
Python = config["tools"]["python"]
Blastp = config["tools"]["Blastp"]
database = config["database"]["bastdb"]
mhchlatype = config["parameter"]["mhchlatype"]
Binsize = config["parameter"]["Binsize"]

print("Binsize",Binsize)


try:
    OutDir = config["OutDir"]
except:
    config["OutDir"] = os.getwd()
    OutDir = config["OutDir"]

rule all:
    input:
        OutDir + "/05.summary/peptide_sequence.csv"

rule EnzymeDigestion:
    output:
        BinStep1 = OutDir + "/01.Enzyme_Digest/bin_step_sequence.fasta",
        BinStep2 = OutDir + "/01.Enzyme_Digest/bin_step_sequence.seq",
        ClipRes = OutDir + "/01.Enzyme_Digest/clip_result.seq",
        DigestPep = OutDir + "/01.Enzyme_Digest/digested_peptide_sequence.fasta"
    shell:
        """
        {Python} {BinCut} --fasta {RawPeptide}  --out1 {output.BinStep1} --out2 {output.BinStep2} --binsize {Binsize}
        {Python} {inhousetoolsconcat} --rule EnzymeDigestion -i {output.BinStep2} -o {output.ClipRes}
        {Python} {GetCutSeq} -f {output.ClipRes} -o {output.DigestPep} --binsize {Binsize}
        """
rule Homology:
    input:
        OutDir + "/01.Enzyme_Digest/digested_peptide_sequence.fasta"
    output:
        blastp =  OutDir + "/02.Homologous_select/peptide_sequence.blastp",
        nonhomo = OutDir + "/02.Homologous_select/non_homologous_peptide_sequence.seq"
    shell:
        """
        {Blastp} -task blastp -db {database} -out {output.blastp} -query {input} -outfmt 6
        {Python} {Homology} -f {output.blastp} -q {input} -i 50 -e 0 -b 100 -o {output.nonhomo} -d {OutDir}/02.Homologous_select -m {mhchlatype} -t {output.blastp}
        """
rule AntigenPresentation:
    input:
        OutDir + "/02.Homologous_select/non_homologous_peptide_sequence.seq"
    output:
        OutDir + "/03.AntigenPresentation/antigen_presentation.seq"
    shell:
        """
        {Python} {inhousetoolsconcat} --rule AntigenPresentation -i {input}  -o {output}
        """
rule Immunogenicity:
    input:
        OutDir + "/02.Homologous_select/non_homologous_peptide_sequence.seq"
    output:
        OutDir + "/04.Immunogenicity/immuno_peptide_sequence.seq"
    shell:
        """
        {Python} {inhousetoolsconcat} --rule Immunogenicity -i {input} -o {output}
        """
rule Summary:
    input:
        antigenpep = OutDir + "/03.AntigenPresentation/antigen_presentation.seq",
        immunopep = OutDir + "/04.Immunogenicity/immuno_peptide_sequence.seq"
    output:
        OutDir + "/05.summary/peptide_sequence.csv"
    shell:
        """
        {Python} {Summary}  -a {input.antigenpep} -i {input.immunopep}  -o {output}
        """

