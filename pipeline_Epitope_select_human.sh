##########################################################################
# File Name:        dev.demo.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Sat 15 Oct 2022 05:05:27 PM CST
##########################################################################
##genelist="MYB=MYB,CDK6=CDK6,IGLL1=IGLL1,RFX8=RFX8,TDT=DNTT,TPOR=MPL,PO4F1=POU4F1,CT451=CT45A1,5HT1F=HTR1F,|WT1=WT1,GPR32,PTX4,F186B=FAM186B"
TAAlist="MYB,CDK6,IGLL1,RFX8,TDT,TPOR,PO4F1,CT451,5HT1F,|WT1,GPR32,PTX4,F186B"
outdir=$1

##########################################################################
##################          config				  ########################
##########################################################################
MHCIlist="HLA-A11:01,HLA-A24:02,HLA-C07:02,HLA-C01:02,HLA-A33:03,HLA-C08:01,\
HLA-C03:04,HLA-A02:01,HLA-B40:01,HLA-C04:01,HLA-B58:01,HLA-B46:01,HLA-B51:01,\
HLA-C03:02,HLA-B38:02,HLA-A02:07,HLA-B15:01,HLA-A02:06,HLA-C03:03,HLA-B15:02,\
HLA-A02:03,HLA-B44:03,HLA-C14:02,HLA-B35:01,HLA-C06:02,HLA-B54:01,HLA-B13:01,\
HLA-B40:02,HLA-B55:02,HLA-A26:01"
MHCIIlist="DRB1_0901,DRB1_1501,DRB1_1202,DRB1_0701,DRB1_0803,DRB1_1101,DRB1_0301,\
DRB1_0405,DRB1_1602,DRB1_1502,DRB1_1302,DRB1_1201,DRB1_1454,DRB1_0406,DRB1_0403,\
DRB1_1405,DRB1_0101,DRB1_1001,DRB1_1301"

#MHCIlist="H-2-Kb,H-2-Db"
#MHCIIlist="H-2-IAb"

netmhcparser=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/parser_for_netmhc.py
bincuter=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/1.1.BinCut.py
blastp=/mnt/data2/wuzengding/03.biotools/software/ncbi-blast-2.13.0+/bin/blastp
proteindb=/mnt/data2/wuzengding/00.database/18.uniprot/index_sp_human_canon_isoform/sp_human_canon_isoform.fasta
homobedmake=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/2.1.HomologySelect_by_outliner_and_makebedfile.py
TMRbedmake=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/5.1.TransMembrance.py
enzymesoft=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/3.1.Concat.py 
homoidcutoff=65

#########################################################################
#################        Pipeline                 #######################
#########################################################################
if  false;then
    echo 'nothing,just used for skipping some steps'
fi
## Stpe1: select protein seuence from uniprot database according selected CTA gene
mkdir -p ${outdir}/01.protein_sequence
echo ${TAAlist}|sed 's:,:\n:g'|while read proteinname;
do  
    echo ${proteinname}
    cat ${proteindb} |\
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
    grep "${proteinname}_"|grep -v "Isoform" |\
    awk -F '\t' '{printf("%s\n%s\n",$1,$2)}' \
    > ${outdir}/01.protein_sequence/protein.${proteinname}.aa.fasta
done

cat ${outdir}/01.protein_sequence/*fasta > \
    ${outdir}/01.protein_sequence/protein.merge.aa.fasta

## Step2:  prediction with netMHCpan for MHC typeI alleles
mkdir -p ${outdir}/02.protein_antigen_prediction
nohup /mnt/data2/wuzengding/03.biotools/software/netMHCpan/netpan41/netMHCpan \
    -a ${MHCIlist} -s \
    -f ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
    -xls \
    -xlsfile ${outdir}/02.protein_antigen_prediction/merge_netMHCpan_out.xls \
    -inptype 0 \
    -BA  \
    > ${outdir}/02.protein_antigen_prediction/merge_netMHCpan_result.xls
nohup /mnt/data2/wuzengding/03.biotools/software/netMHCpan/netpanii40/netMHCIIpan \
	-a ${MHCIIlist} -s \
	-f ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
	-xls \
	-xlsfile ${outdir}/02.protein_antigen_prediction/merge_netMHCIIpan_out.xls \
	-inptype 0 \
	-BA \
	> ${outdir}/02.protein_antigen_prediction/merge_netMHCIIpan_result.xls


python ${netmhcparser} \
	-u ${outdir}/02.protein_antigen_prediction/merge_netMHCIIpan_result.xls \
	-i ${outdir}/02.protein_antigen_prediction/merge_netMHCpan_result.xls \
    -o ${outdir}/02.protein_antigen_prediction \
    -b Y

## Step3: prediction of homology
mkdir -p ${outdir}/03.homologous
python3 ${bincuter} \
    -f ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
    -b 12  \
    -u ${outdir}/03.homologous/merge_bin_step_sequence.seq  \
    -o ${outdir}/03.homologous/merge_bin_step_sequence.fasta
 
${blastp} -task blastp \
    -db ${proteindb} \
    -out ${outdir}/03.homologous/peptide_sequence.blastp \
    -query ${outdir}/03.homologous/merge_bin_step_sequence.fasta \
    -outfmt 6

python3 ${homobedmake} \
    -f ${outdir}/03.homologous/peptide_sequence.blastp \
    -i ${homoidcutoff} \
    -o ${outdir}/03.homologous

## make blastp with fasta generated from netMHC.csv 
tail -n+2 ${outdir}/02.protein_antigen_prediction/parsed_res_MHCI/merge.netMHCI.csv |\
    awk -F',' '{print ">"$11"_"$1","$3}'|sort|uniq|tr ',' '\n' > \
        ${outdir}/03.homologous/protein.netMHC.aa.fasta
tail -n+2 ${outdir}/02.protein_antigen_prediction/parsed_res_MHCII/merge.netMHCII.csv|\
    awk -F',' '{print ">"$7"_"$1","$3}'|sort|uniq|tr ',' '\n' >> \
        ${outdir}/03.homologous/protein.netMHC.aa.fasta
${blastp} -task blastp \
    -db ${proteindb} \
    -out ${outdir}/03.homologous/peptide_netMHC.blastp \
    -query ${outdir}/03.homologous/protein.netMHC.aa.fasta \
    -outfmt 6
    
## Step4: transmembrane prediction
cd  ${outdir}
biolib run DTU/DeepTMHMM --fasta ${outdir}/01.protein_sequence/protein.merge.aa.fasta
mv biolib_results  04.TransMembrane.DeepTMHMM


python3  ${TMRbedmake} \
    -f ${outdir}/04.TransMembrane.DeepTMHMM/TMRs.gff3 \
    -o ${outdir}/04.TransMembrane.DeepTMHMM

## Step5: enzyme digestions prediction
mkdir -p ${outdir}/05.EnzymeDigest
python enzymesoft 
    -r EnzymeDigestion \
    -i ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
    -o ${outdir}/05.EnzymeDigest/protein.merge.enzymedigest

## Step5: Summary and deliverables
mkdir -p ${outdir}/06.Deliverables 
    cp ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
                ${outdir}/06.Deliverables/01.Protein.merge.aa.fasta
    cp ${outdir}/02.protein_antigen_prediction/parsed_res_MHCI/Immunogenicity.netMHCI.bed \
                ${outdir}/06.Deliverables/02.Immunogenicity.netMHCI.bed
    cp ${outdir}/02.protein_antigen_prediction/parsed_res_MHCII/Immunogenicity.netMHCII.bed \
                ${outdir}/06.Deliverables/02.Immunogenicity.netMHCII.bed
    cp ${outdir}/03.homologous/homo_peptide.bed  \
                ${outdir}/06.Deliverables/03.Homo_peptide.bed
    cp ${outdir}/04.TransMembrane.DeepTMHMM/TransMembrane.bed \
                ${outdir}/06.Deliverables/04.TransMembrane.bed

echo ${TAAlist}|sed 's:,:\n:g'|while read proteinname;
do 
    echo ${proteinname}
    proteinname_upp=$(echo ${proteinname}|tr  [:lower:] [:upper:])
    epitopname=$(grep ${proteinname_upp} \
        ${outdir}/06.Deliverables/02.Immunogenicity.netMHCI.bed |\
    cut -f1 |uniq)
    
    echo ${epitopname}
    ## trans seqid of file of '01.Protein.merge.aa.fasta' 
    seqid=$(grep ${proteinname} ${outdir}/06.Deliverables/01.Protein.merge.aa.fasta)
    sed -i "s;${seqid};\>${epitopname};g" ${outdir}/06.Deliverables/01.Protein.merge.aa.fasta
    
    ## trans seqid of file of '02.Immunogenicity.netMHCII.bed' 
    echo "trans seqid of file of '02.Immunogenicity.netMHCII.bed'"
    grep ${proteinname_upp} ${outdir}/06.Deliverables/02.Immunogenicity.netMHCII.bed |\
    cut -f1 |uniq|while read seqid;
    do  
        echo ${seqid}
        sed -i "s;${seqid};${epitopname};g" ${outdir}/06.Deliverables/02.Immunogenicity.netMHCII.bed
    done
    
    ## trans seqid of file of '03.Homo_peptide.bed' 
    echo "trans seqid of file of 'homo_peptide.bed"
    grep ${proteinname} ${outdir}/06.Deliverables/03.Homo_peptide.bed |\
    cut -f1 |uniq| while read seqid;
    do  
        echo ${seqid}|tr '\|' '_'
        sed -i "s;${seqid};${epitopname};g" ${outdir}/06.Deliverables/03.Homo_peptide.bed
    done
    
    ## trans seqid of file of '04.TransMembrane.bed' 
    echo "trans seqid of file of '04.TransMembrane.bed"
    grep ${proteinname_upp} ${outdir}/06.Deliverables/04.TransMembrane.bed |\
    cut -f1 |uniq| while read seqid;
    do  
        echo ${seqid}
        sed -i "s;${seqid};${epitopname};g" ${outdir}/06.Deliverables/04.TransMembrane.bed
    done
done
