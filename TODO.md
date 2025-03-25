


# preprocessign only
kraken-tools preprocess \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --output-dir results/preprocessed/ \
    --threads 8


# taxonomic classification only
kraken-tools classify \
    --input-fastq clean_reads.fastq \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --output-dir results/taxonomy/ \
    --taxonomic-level S \
    --threads 8

# process existing kraken or bracken files
kraken-tools process \
    --kreport-dir kraken_reports/ \
    --bracken-dir bracken_files/ \
    --sample-key metadata.csv \
    --output-dir results/processed/ \
    --min-abundance 0.01 \
    --min-prevalence 0.1


# Downstream analysis
kraken-tools analyze \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/analysis/ \
    --group-col "Group"

# differential abundance testing
kraken-tools diff-abundance \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/diff_abundance/ \
    --group-col "Group" \
    --methods aldex2,ancom,ancom-bc

# GLMM analysis
kraken-tools glmm \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/glmm/ \
    --formula "Count ~ Group + (1|Subject)" \
    --model negbin





# all steps
kraken-tools full-pipeline \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "Group" \
    --min-abundance 0.01 \
    --min-prevalence 0.1 \
    --threads 8