##########################################################################
# File Name:     pipeline_Epitope_analysis.sh
# Author:		 wuzengding
# mail:			 wuzengding@163.com
# Created Time:	 Sat 15 Oct 2022 05:05:27 PM CST
# Revise History:	 Sat 17 Mar 2023 05:05:27 PM CST-cuifenglei
# Revise History:     17 Jun 2024   wuzengding
##########################################################################
#! /usr/bin/bash

#!/bin/bash

# 打印使用帮助信息
usage() {
  echo "Usage: $0 
      -f <faa_file> -l <HLA_genotypes> -o <output_dir> -p <output_prefix>"
  exit 1
}


# 解析参数
while getopts "f:l:o:p:" opt; do
  case ${opt} in
    f )
      faa_file=$OPTARG
      ;;
    l )
      HLA_genotypes=$OPTARG
      ;;
    o )
      output_dir=$OPTARG
      ;;
    p )
      output_prefix=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done


# 检查文件是否已提供
if [ -z "${faa_file}" ] || [ -z "${HLA_genotypes}" ] || [ -z "${output_dir}" ] || [ -z "${output_prefix}" ]; then
  usage
fi

output_dir=$(realpath ${output_dir})
# 打印解析后的参数
echo "faa_file: ${faa_file}"
echo "HLA_genotypes: ${HLA_genotypes}"
echo "output_dir: ${output_dir}"
echo "output_prefix: ${output_prefix}"


# 示例：使用解析后的参数执行命令
# command_to_run -f "${faa_file}"  -h "${HLA_genotypes}" -o "${output_dir}/${output_prefix}_result"


##########################################################################
##################		  config				  ########################
##########################################################################
netMHCpan=/mnt/user/wzd/03.biotools/software/netMHCpan/netpan41/netMHCpan
netMHCIIpan=/mnt/user/wzd/03.biotools/software/netMHCpan/netpanii40/netMHCIIpan
blastp=/mnt/user/wzd/03.biotools/software/ncbi-blast-2.13.0+/bin/blastp
basedir=$(cd "$(dirname "$0")" && pwd)


####################       tools                #######################
# 定义运行NetMHCpan函数
run_netMHCpan() {
  local faa_file=$1
  local mhci_genotypes=$2
  local output_dir=$3
  local output_prefix=$4

  ${netMHCpan} \
    -a ${mhci_genotypes} -s \
    -f ${faa_file} \
    -xls \
    -xlsfile ${output_dir}/${output_prefix}_netMHCpan_out.xls \
    -inptype 0 \
    -BA \
    > ${output_dir}/${output_prefix}_netMHCpan_result.xls
}

# 定义运行NetMHCIIpan函数
run_netMHCIIpan() {
  local faa_file=$1
  local mhcii_genotypes=$2
  local output_dir=$3
  local output_prefix=$4

  ${netMHCIIpan} \
    -a ${mhcii_genotypes} -s \
    -f ${faa_file} \
    -xls \
    -xlsfile ${output_dir}/${output_prefix}_netMHCIIpan_out.xls \
    -inptype 0 \
    -BA \
    > ${output_dir}/${output_prefix}_netMHCIIpan_result.xls
}


##################   step0: check MHC genetype,species  #########################
json_output=$(python3 ${basedir}/script/identify_mhc_genotype.py  \
            -g "$HLA_genotypes" \
            --mhci-file  ${basedir}/lib/MHCI_genetype_list \
            --mhcii-file ${basedir}/lib/MHCII_genetype_list)



# 使用jq解析JSON输出
mhci_genotypes=$(echo $json_output | jq -r '.mhc_genotypes."MHC I" | join(",")')
mhcii_genotypes=$(echo $json_output | jq -r '.mhc_genotypes."MHC II" | join(",")')

# 解析species_genotypes
species_genotypes=$(echo $json_output | jq -r '.species_genotypes | to_entries | map("\(.key): \(.value | join(", "))") | join("; ")' )
species=$(echo $json_output | jq -r '.species_genotypes | to_entries | map("\(.key)")| join("; ")')
  ## species_genotypes eg: human: DRB1_1412, HLA-A01:127; mouse: H-2-Kb
  ## species eg:human;mouse
species_count=$(echo $species | tr ';' '\n' | wc -l)
if [ "$species_count" -ge 2 ]; then
    echo "Warning: More than one species found: $species"
    exit 1
fi

# 输出解析结果
echo "MHC I genotypes: $mhci_genotypes"
echo "MHC II genotypes: $mhcii_genotypes"
echo -e "Species genotypes:\n$species_genotypes"

echo "$(realpath ${output_dir}/01.protein_sequence)"





##################  step1: Pre-precess fasta file    ##########################

if [ -e "${output_dir}/01.protein_sequence" ]; then
    echo "Step1: protein sequence  existed, the fasta id format translating  step passed!"
else
    echo "Step1: fasta id format translating....."
    mkdir -p "${output_dir}/01.protein_sequence"
    #cp ${faa_file}  ${output_dir}/01.protein_sequence
    
    python3 ${basedir}/script/preprocess_fasta.py \
      -f ${faa_file} \
      -o ${output_dir}/01.protein_sequence \
      -p ${output_prefix}
fi

processed_var_faa="${output_dir}/01.protein_sequence/${output_prefix}_var_seq.faa"
processed_ref_faa="${output_dir}/01.protein_sequence/${output_prefix}_ref_seq.faa"


######### Step2:  prediction with netMHCpan for MHC typeI/II alleles ###########

run_netmhc() {
    local sequence_file=$1
    local output_subdir=$2
    local output_prefix=$3
    local mhci_genotypes=$4
    local mhcii_genotypes=$5
    local bindlevel=$6

    echo "Running netMHC predictions for ${sequence_file} ..."

    mkdir -p "${output_dir}/${output_subdir}"
    mkdir -p "${output_dir}/${output_subdir}/parsed"

    if [ -n "${mhci_genotypes}" ]; then
        run_netMHCpan "${sequence_file}" "${mhci_genotypes}" "${output_dir}/${output_subdir}" "${output_prefix}"
        python3 ${basedir}/script/parser_netmhc.py \
            -t "netMHCpan" \
            -i "${output_dir}/${output_subdir}/${output_prefix}_netMHCpan_result.xls" \
            -o "${output_dir}/${output_subdir}/parsed" \
            -p "${output_prefix}" \
            -y "csvfasta" \
            -b "${bindlevel}"
            
    fi

    if [ -n "${mhcii_genotypes}" ]; then
        run_netMHCIIpan "${sequence_file}" "${mhcii_genotypes}" "${output_dir}/${output_subdir}" "${output_prefix}"
        python3 ${basedir}/script/parser_netmhc.py \
            -t "netMHCIIpan" \
            -i "${output_dir}/${output_subdir}/${output_prefix}_netMHCIIpan_result.xls" \
            -o "${output_dir}/${output_subdir}/parsed" \
            -p "${output_prefix}" \
            -y "csvfasta" \
            -b "${bindlevel}"
    fi
}

if [ -e "${output_dir}/02.protein_antigen_prediction_var" ]; then
    echo "Step2: protein_antigen_prediction var existed, this step passed!"
else
    if [ -e ${processed_var_faa} ];then
        echo "Step2: prediction with netMHC for var seq Starting......"
        run_netmhc "${processed_var_faa}" "02.protein_antigen_prediction_var" "${output_prefix}" "${mhci_genotypes}" "${mhcii_genotypes}" "SBWB"
    else
        echo " ${processed_var_faa} non-exist!"
        exit 1
    fi
fi

if [ -e "${output_dir}/02.protein_antigen_prediction_ref" ]; then
    echo "Step2.2: protein_antigen_prediction ref existed, this step passed!"
else
    if [ -e ${processed_ref_faa} ];then
        echo "Step2.2: prediction with netMHC for ref seq Starting......"
        run_netmhc "${processed_ref_faa}" "02.protein_antigen_prediction_ref" "${output_prefix}" "${mhci_genotypes}" "${mhcii_genotypes}" "all"
    else
        echo " ${processed_ref_faa} non-exist!"
        exit 1
    fi
fi


############      Step3: prediction of homology        ###################
if [ -e "${output_dir}/03.homologous" ]; then
    echo "Step3: 03.homologous existed, this step passed!"
else
    echo "Step3: homologous analysis Starting......"
    mkdir -p "${output_dir}/03.homologous"
    
    if [ ${species} = "mouse" ];then
        proteindb="${basedir}/lib/index_sp_mouse_blastdb/sp_mouse_canon_isoform.fasta"
    elif [ ${species} = "human" ];then
        proteindb="${basedir}/lib/index_sp_human_blastdb/sp_human_canon_isoform.fasta"
    else
        echo "Warning: your species is not mouse or human"
        exit 1
    fi
    
    if [ -n "${mhci_genotypes}" ]; then
        ${blastp} -task blastp \
            -db ${proteindb} \
            -out ${output_dir}/03.homologous/${output_prefix}_netMHCpan.blastp \
            -query "${output_dir}/02.protein_antigen_prediction_var/parsed/${output_prefix}_netMHCpan.faa" \
            -outfmt 6

        python3 ${basedir}/script/parser_blastp.py \
            -b ${output_dir}/03.homologous/${output_prefix}_netMHCpan.blastp \
            -f ${proteindb} \
            -o ${output_dir}/03.homologous/ \
            -p ${output_prefix}_netMHCpan \
            -c "0.8"  ## pident cutoff

        homo_faa="${output_dir}/03.homologous/${output_prefix}_netMHCpan_homologous.faa"
        run_netMHCpan "${homo_faa}" "${mhci_genotypes}" "${output_dir}/03.homologous" "${output_prefix}_homologous"
        
        python3 ${basedir}/script/parser_netmhc.py \
            -t "netMHCpan" \
            -i "${output_dir}/03.homologous/${output_prefix}_homologous_netMHCpan_result.xls" \
            -o "${output_dir}/03.homologous/parsed" \
            -p "${output_prefix}_homologous" \
            -y "csv"     
    fi
        
    if [ -n "${mhcii_genotypes}" ]; then
        ${blastp} -task blastp \
            -db ${proteindb} \
            -out ${output_dir}/03.homologous/${output_prefix}_netMHCIIpan.blastp \
            -query "${output_dir}/02.protein_antigen_prediction_var/parsed/${output_prefix}_netMHCIIpan.faa" \
            -outfmt 6

        python3 ${basedir}/script/parser_blastp.py \
            -b ${output_dir}/03.homologous/${output_prefix}_netMHCIIpan.blastp \
            -f ${proteindb} \
            -o ${output_dir}/03.homologous/ \
            -p ${output_prefix}_netMHCIIpan \
            -c "0.8"  ## pident cutoff
            
        homo_faa="${output_dir}/03.homologous/${output_prefix}_netMHCIIpan_homologous.faa"    
        run_netMHCIIpan "${homo_faa}" "${mhcii_genotypes}" "${output_dir}/03.homologous" "${output_prefix}_homologous"
        
        python3 ${basedir}/script/parser_netmhc.py \
            -t "netMHCIIpan" \
            -i "${output_dir}/03.homologous/${output_prefix}_homologous_netMHCIIpan_result.xls" \
            -o "${output_dir}/03.homologous/parsed" \
            -p "${output_prefix}_homologous" \
            -y "csv"   
    fi

    
fi


###########       Step4: transmembrane prediction       ##############
processed_var_faa=$(realpath ${processed_var_faa})
if [ -e ${output_dir}/04.TransMembrane.DeepTMHMM ];then
    echo "Step4: TransMembrane.DeepTMHMM existed, this step passed!"
else
    echo "Step4: TransMembrane.DeepTMHMM  Starting......"
    #echo "biolib run DTU/DeepTMHMM --fasta ${outdir}/01.protein_sequence/protein.merge.aa.fasta"
    cd  ${output_dir}
    biolib run DTU/DeepTMHMM --fasta  ${processed_var_faa} < /dev/null
    if [ -d  "$output_dir/04.TransMembrane.DeepTMHMM/biolib_results" ];then
        rm -rf $output_dir/04.TransMembrane.DeepTMHMM/biolib_results
    fi
    mv biolib_results  04.TransMembrane.DeepTMHMM
        
    python3  ${TMRbedmake} \
        -f ${output_dir}/04.TransMembrane.DeepTMHMM/TMRs.gff3 \
        -o ${output_dir}/04.TransMembrane.DeepTMHMM
fi


###########      Step5: enzyme digestions prediction   ###############

echo ${output_dir}
cd ${output_dir}
if [ -e ${output_dir}/05.EnzymeDigest ];then
    echo "Step5:  EnzymeDigest existed, this step passed!"
else
    echo "Step5:  EnzymeDigest Starting..........."
    mkdir -p ${output_dir}/05.EnzymeDigest
    true >${output_dir}/05.EnzymeDigest/protein.merge.enzymedigest
    #python ${enzymesoft} \
    #	-r EnzymeDigestion \
    #	-i ${outdir}/01.protein_sequence/protein.merge.aa.fasta \
    #	-o ${outdir}/05.EnzymeDigest/protein.merge.enzymedigest
fi

###############     summary results             ##################
echo $(dirname ${output_dir})
python ${basedir}/script/summarize.py \
    -d ${output_dir} \
    -o ${output_dir}/Deliverable \
    -p ${output_prefix} \
    -n ${output_prefix} \
    -m 'all'