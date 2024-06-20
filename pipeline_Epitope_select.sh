##########################################################################
# File Name:        dev.demo.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Sat 15 Oct 2022 05:05:27 PM CST
##########################################################################
TAAlist="CT83,MAGEA3"
outdir=$1

##########################################################################
##################          config				  ########################
##########################################################################
MHCIlist="HLA-A11:01,HLA-A24:02,HLA-C07:02,HLA-C01:02,HLA-A33:03,HLA-C08:01,\
HLA-C03:04,HLA-A02:01,HLA-B40:01,HLA-C04:01,HLA-B58:01,HLA-B46:01,HLA-B51:01,\
HLA-C03:02,HLA-B38:02,HLA-A02:07,HLA-B15:01,HLA-A02:06,HLA-C03:03,HLA-B15:02,\
HLA-A02:03,HLA-B44:03,HLA-C14:02,HLA-B35:01,HLA-C06:02,HLA-B54:01,HLA-B13:01,\
HLA-B40:02,HLA-B55:02,HLA-A26:01"

netmhcparser=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/parser_for_netmhc.py
bincuter=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/1.1.BinCut.py
blastp=/mnt/data2/wuzengding/03.biotools/software/ncbi-blast-2.13.0+/bin/blastp
proteindb=/mnt/data2/wuzengding/00.database/18.uniprot/index_sp_human_canon_isoform/sp_human_canon_isoform.fasta
homobedmake=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/2.1.HomologySelect_by_outliner_and_makebedfile.py
TMRbedmake=/mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/pipeline/5.1.TransMembrance.py
homoidcutoff=65

#########################################################################
#################        Pipeline                 #######################
#########################################################################
if  false;then
    echo 'nothing'
## Stpe1: select protein seuence from uniprot database according selected CTA gene
mkdir -p ${outdir}/01.protein_sequence
echo ${TAAlist}|sed 's:,:\n:g'|while read proteinname;
do  
    echo ${proteinname}
    cat ${proteindb} |\
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |grep ${proteinname}|\
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

python ${netmhcparser} -i \
    ${outdir}/02.protein_antigen_prediction/merge_netMHCpan_result.xls \
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

## Step4: transmembrane prediction
cd  ${outdir}
biolib run DTU/DeepTMHMM --fasta ${outdir}/01.protein_sequence/protein.merge.aa.fasta
mv biolib_results  04.TransMembrane.DeepTMHMM

fi
python3  ${TMRbedmake} \
    -f ${outdir}/04.TransMembrane.DeepTMHMM/TMRs.gff3 \
    -o ${outdir}/04.TransMembrane.DeepTMHMM

## Step5: Summary and deliverables
mkdir -p ${outdir}/05.Deliverables 
    cp ${outdir}/01.protein_sequence/protein.merge.aa.fasta ${outdir}/05.Deliverables/
    cp ${outdir}/02.protein_antigen_prediction/*/*.bed ${outdir}/05.Deliverables/
    cp ${outdir}/03.homologous/homo_peptide.bed ${outdir}/05.Deliverables/
    cp ${outdir}/04.TransMembrane.DeepTMHMM/TransMembrane.bed ${outdir}/05.Deliverables/