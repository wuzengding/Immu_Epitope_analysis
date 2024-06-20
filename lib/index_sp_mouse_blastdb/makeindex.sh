##########################################################################
# File Name:        makeindex.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 09 May 2022 04:16:29 PM CST
##########################################################################

/mnt/data2/wuzengding/03.biotools/software/ncbi-blast-2.13.0+/bin/makeblastdb -in /mnt/data2/wuzengding/05.pipeline_dev/Antigen_Design/01.mouse_protein_sequence_data/sp_mouse_canon_isoform.fasta -parse_seqids -blastdb_version 5 -title 'mouse_protein' -dbtype prot -out /mnt/data2/wuzengding/00.database/index_mus_blastdb/index_sp_mouse_canon_isoform_blastdb/sp_mouse_canon_isoform.fasta
