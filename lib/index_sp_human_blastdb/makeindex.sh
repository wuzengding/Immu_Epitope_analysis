##########################################################################
# File Name:        makeindex.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 17 Oct 2022 05:01:28 PM CST
##########################################################################

/mnt/data2/wuzengding/03.biotools/software/ncbi-blast-2.13.0+/bin/makeblastdb  -in /mnt/data2/wuzengding/00.database/18.uniprot/index_sp_human_canon_isoform/sp_human_canon_isoform.fasta  -parse_seqids -blastdb_version 5 -title 'human_protein' -dbtype prot
