# **how to extract out expression/reads from fastq files**

grab fasta files files or convert genes in genbank files to fasta python parsse_gb_2_lwgv.py --\> generates {name}\_gene.fasta

generate index via blast or kallisto## bash mk_blast_idx.sh --\> .idx/{name} bash mk_kall_idx.sh --\> .idx/{name}

map reads to index\$\$ bash run_kallisto.sh --\> {organism}/abundance.tsv, {organism}/abundance.h5, {organism}/run_info.json

extract expressions and generate files --\> exp/{organism}\_table.csv schema shows how to organize expressions are extracted out based on organism/group and pattern
