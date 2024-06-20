#!/bin/bash

# 定义一个关联数组来映射MHC基因型前缀到对应的物种名称
declare -A prefix_to_species=(
    ["BoLA"]="bovine"
    ["Chi"]="chicken"
    ["DLA"]="dog"
    ["Eqca"]="equine"
    ["Gogo"]="gorilla"
    ["H-2"]="mouse"
    ["H2"]="mouse"
    ["HLA"]="human"
    ["DRB1"]="human"
    ["DRB3"]="human"
    ["DRB4"]="human"
    ["DRB5"]="human"
    ["Mamu"]="macaque"
    ["Patr"]="chimpanzee"
    ["SLA"]="swine"
)

# 处理 netMHCpan -listMHC 的输出
/mnt/user/wzd/03.biotools/software/netMHCpan/netpan41/netMHCpan -listMHC | grep -vE '^#|^$' | while IFS= read -r line; do
    matched=false
    for prefix in "${!prefix_to_species[@]}"; do
        if [[ "$line" == "$prefix"* ]]; then
            echo "${prefix_to_species[$prefix]}|$line"
            matched=true
            break
        fi
    done
    if [[ "$matched" != true ]]; then
        echo "unknown|$line"
    fi
done > ./MHCI_genetype_list

# 处理 netMHCIIpan -list 的输出
/mnt/user/wzd/03.biotools/software/netMHCpan/netpanii40/netMHCIIpan -list | grep -vE '^#|^$' | while IFS= read -r line; do
    matched=false
    for prefix in "${!prefix_to_species[@]}"; do
        if [[ "$line" == "$prefix"* ]]; then
            echo "${prefix_to_species[$prefix]}|$line"
            matched=true
            break
        fi
    done
    if [[ "$matched" != true ]]; then
        echo "unknown|$line"
    fi
done > ./MHCII_genetype_list