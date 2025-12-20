pdb2fasta.py

mmseqs createdb all_sequences.fasta mmseqs_DB

mmseqs cluster mmseq_db/mmseqs_DB mmseqs_DB_clu tmp_dir --min-seq-id 0.3 -c 0.8 --cov-mode 0

mmseqs createtsv mmseq_db/mmseqs_DB mmseq_db/mmseqs_DB mmseqs_DB_clu cluster_results.tsv

awk '{print $1}' cluster_results.tsv | sort | uniq > representative_pdbs.txt

cut -d '_' -f 1 representative_pdbs.txt | sort | uniq > representative_pdbs_clean.txt
