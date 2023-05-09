#!/usr/bin/env python3

from pathlib import Path


#############
# FUNCTIONS #
#############

def get_guppy_fastq_files(wildcards):
    # need to handle old version of guppy that don't have pass/fail dirs
    if Path('output/010_basecall/guppy/pass').is_dir():
        return('output/010_basecall/guppy/pass/{read}.fastq')
    else:
        return('output/010_basecall/guppy/{read}.fastq')


def aggregate_reads(wildcards):
    idlist = checkpoints.generate_read_id_list.get(**wildcards).output['idlist']
    with open(idlist, 'rt') as f:
        read_ids = [line.rstrip() for line in f]
    return(
        snakemake.io.expand(
            'output/020_porechop/{read}.fastq',
            read=read_ids))


###########
# GLOBALS #
###########

# should be a symlink to the fast5 directory
# TODO: handle zip files
fast5_path = Path('data/reads').resolve()

# CONTAINERS
bioconductor = 'docker://ghcr.io/tomharrop/r-containers:bioconductor_3.17'
biopython = 'docker://quay.io/biocontainers/biopython:1.78'
guppy = 'docker://ghcr.io/tomharrop/container-guppy:v6.4.6_cv2'
minionqc = 'shub://TomHarrop/ont-containers:minionqc_1.4.1'
pigz = 'docker://quay.io/biocontainers/pigz:2.3.4'
porechop = 'docker://quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'


#########
# RULES #
#########

rule target:
    input:
        'output/015_qc/guppy/summary.yaml',
        'output/015_qc/read_distribution.pdf',
        'output/processed_reads.fastq.gz'

# process reads
rule compress_processed_reads:
    input:
        'output/020_porechop/processed_reads.fastq'
    output:
        'output/processed_reads.fastq.gz'
    log:
        'output/logs/compress_processed_reads.log'
    threads:
        10
    resources:
        time = 20
    container:
        pigz
    shell:
        'pigz -p {threads} '
        '-9 '
        '<{input} '
        '>{output} '
        '2> {log}'


rule aggregate_reads:
    input:
        aggregate_reads
    output:
        pipe('output/020_porechop/processed_reads.fastq')
    shell:
        'cat {input} > {output}'

rule porechop:
    input:
        'output/010_basecall/tmp/{read}.fastq'
    output:
        temp('output/020_porechop/{read}.fastq')
    log:
        'output/logs/porechop/{read}.log'
    threads:
        1
    resources:
        time = 10
    container:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--discard_middle '
        '&> {log}'


# extract the guppy output for processing
rule unzip_fastq_file:
    input:
        'output/010_basecall/compressed_reads/{read}.fastq.gz'
    output:
        temp('output/010_basecall/tmp/{read}.fastq')
    log:
        'output/logs/unzip_fastq_file/{read}.log'
    container:
        pigz
    shell:
        'pigz -d <{input} >{output} 2>{log}'

# compress the guppy output for storage
rule gzip_fastq_file:
    input:
        read = get_guppy_fastq_files,
        #  make sure basecalling has happened
        summary_file = 'output/010_basecall/guppy/sequencing_summary.txt'
    output:
        'output/010_basecall/compressed_reads/{read}.fastq.gz'
    log:
        'output/logs/gzip_fastq_file/{read}.log'
    threads:
        10
    resources:
        time = 1
    container:
        pigz
    shell:
        'pigz -p {threads} '
        '-9 '
        '<{input.read} '
        '>{output} '
        '&& rm {input.read} '
        '&> {log}'


# qc
rule plot_read_distribution:
    input:
        seqsum = 'output/010_basecall/guppy/sequencing_summary.txt'
    output:
        plot = 'output/015_qc/read_distribution.pdf'
    threads:
        1
    container:
        bioconductor
    log:
        'output/logs/plot_read_distribution.log'
    script:
        'src/plot_read_distribution.R'

rule minionqc:
    input:
        'output/010_basecall/guppy/sequencing_summary.txt'
    output:
        'output/015_qc/guppy/summary.yaml'
    params:
        search_dir = 'output/010_basecall/guppy',
        outdir = 'output/015_qc'
    threads:
        1
    resources:
        mem_mb = 8000
    container:
        minionqc
    log:
        'output/logs/minionqc.log'
    shell:
        'MinIONQC.R '
        '--processors={threads} '
        '--input={params.search_dir} '
        '--outputdirectory={params.outdir} '
        '&> {log}'


# generate a list of read IDs
checkpoint generate_read_id_list:
    input:
        seqsum = 'output/010_basecall/guppy/sequencing_summary.txt'
    output:
        idlist = 'output/010_basecall/read_list.txt'
    threads:
        1
    resources:
        time = 1
    log:
        'output/logs/generate_read_id_list.log'
    container:
        biopython
    script:
        'src/generate_read_id_list.py'

rule guppy:
    input:
        fast5_path
    output:
        'output/010_basecall/guppy/sequencing_summary.txt',
        # Marking the fail directory as temp means Snakemake will delete it
        # after basecalling completes. This is to save storage space.
        f = temp(directory(f'output/010_basecall/guppy/fail'))
    params:
        outdir = 'output/010_basecall/guppy',
        config = 'dna_r9.4.1_450bps_sup.cfg'
    log:
        'output/logs/guppy.log'
    threads:
        3
    resources:
        partition = 'gpu-a100',
        gres_flag = '--gres=gpu:1',
        time = 480 * 5,
        mem_mb = 40000
    container:
        guppy
    shell:
        'guppy_basecaller '
        '--device auto '        # enable GPU
        '--input_path {input} '
        '--save_path {params.outdir} '
        '--config {params.config} '
        '--verbose_logs '
        '--recursive '
        '&> {log}'

