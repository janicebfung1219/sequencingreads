gunzip -c fa/mutans_genes.fa.gz  |  makeblastdb  -input_type fasta -dbtype nucl -parse_seqids  -in - -title mutans -out idx/mutans&


