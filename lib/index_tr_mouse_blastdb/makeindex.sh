##########################################################################
# File Name:        makeindex.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 09 May 2022 04:16:29 PM CST
##########################################################################

/mnt/data2/wuzengding/03.biotools/software/ncbi-blast-2.13.0+/bin/makeblastdb -in /mnt/data2/wuzengding/00.database/18.uniprot/tr_mouse_canon_isoform.fasta -out /mnt/data2/wuzengding/00.database/index_mus_blastdb/index_tr_mouse_canon_isoform_blastdb/tr_mouse_canon_isoform.fasta  -parse_seqids -blastdb_version 5 -title 'mouse_protein' -dbtype prot
