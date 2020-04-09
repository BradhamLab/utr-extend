from glob import glob
import os
import utils

configfile: 'config.yaml'
shell.prefix('source activate alignment; ')

lib_match = config['params']['library_regex']
# change this
LIBRARIES = [os.path.basename(x).split('.')[0]\
             for x in glob(os.path.join(config['dir']['data'] + f'*{lib_match}*'))]
# LIBRARIES = ['L001', 'L002', 'L003', 'L004']

rule all:
    input:
        os.path.join(config['dir']['out'], 'results',
                    'filtered.matched.bedtools.closest.out'),
        os.path.join(config['dir']['out'], 'results',
                    'unmatched.bedtools.closest.out'),

rule gtf_to_bed:
    input:
        gtf=config['genome']['gtf']
    output:
        bed=temp(os.path.join(config['dir']['out'], 'annotations',
                              'annotations.bed'))
    shell:
        'gtf2bed < {input.gtf} > {output.bed}'

rule exon_annos:
    input:
        bed=os.path.join(config['dir']['out'], 'annotations', 'annotations.bed')
    output:
        bed=temp(os.path.join(config['dir']['out'], 'annotations',
                'exon.annotations.bed'))
    shell:
        "awk '{{if ($8 ~ /exon/) {{ print }} }}' {input} > {output}"

rule sort_annotations:
    input:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'exon.annotations.bed')
    output:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'sorted.exon.annotations.bed')
    shell:
        'sort -k1,1 -k2,2n {input.bed} > {output.bed}'

rule build_star_index:
    input:
        fasta=config['genome']['fasta'],
        gtf=config['genome']['gtf'],
    params:
        index_dir=config['STAR']['index'],
        chr_n_bits=utils.estimate_STAR_ChrBinNbits(config['genome']['fasta'], 60),
    output:
        os.path.join(config['STAR']['index'], 'Genome')
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 60 --genomeChrBinNbits {params.chr_n_bits}"

rule combine_fastq:
    input:
        fastq=expand(os.path.join(config['dir']['data'], '{library}' + '.fastq'),
                     library=LIBRARIES)
    output:
        fastq=os.path.join(config['dir']['out'], 'fastq', 'combined.fastq')
    shell:
        "cat {input.fastq} > {output.fastq}"

rule align_3_prime_utr_reads:
    input:
        index=os.path.join(config['STAR']['index'], 'Genome'),
        fastq=os.path.join(config['dir']['out'], 'fastq', 'combined.fastq')
    output:
        sam=temp(os.path.join(config['dir']['out'], "STAR", "Aligned.out.sam"))
    params:
        index=os.path.join(config['STAR']['index']),
        prefix=os.path.join(config['dir']['out'], "STAR") + '/'
    shell:
        'STAR --runMode alignReads --outSAMtype SAM --readFilesCommand zcat '
        '--genomeDir {params.index} --outFileNamePrefix {params.prefix} '
        '--readFilesIn {input.fastq} --outReadsUnmapped FASTX'

rule sam_to_bed:
    input:
        sam=os.path.join(config['dir']['out'], "STAR", "Aligned.out.sam")
    output:
        bed=temp(os.path.join(config['dir']['out'], 'bed', 'aligned.bed'))
    shell:
        'sam2bed < {input.sam} > {output.bed}'

rule sort_alignment_bed:
    input:
        bed=os.path.join(config['dir']['out'], 'bed', 'aligned.bed')
    output:
        bed=os.path.join(config['dir']['out'], 'bed', 'aligned.sorted.bed')
    shell:
        'sort -k1,1 -k2,2n {input.bed} > {output.bed}'

rule remove_matched:
    input:
        bed=os.path.join(config['dir']['out'], 'bed', 'aligned.sorted.bed'),
        ref=os.path.join(config['dir']['out'], 'annotations',
                         'sorted.exon.annotations.bed')
    output:
        bed=os.path.join(config['dir']['out'], 'bed', 'unannotated.bed')
    shell:
        'bedtools subtract -a {input.bed} -b {input.ref} -s -A > {output.bed}'

rule merge_missed_alignments:
    input:
        bed=os.path.join(config['dir']['out'], 'bed', 'unannotated.bed')
    output:
        bed=temp(os.path.join(config['dir']['out'], 'bed', 'merged.bed'))
    shell:
        'bedtools merge -i {input.bed} -s -c 4,5,6 -o count,distinct,distinct'

rule find_closest:
    input:
        ref=os.path.join(config['dir']['out'], 'annotations',
                         'sorted.exon.annotations.bed'),
        bed=os.path.join(config['dir']['out'], 'bed',
                         'merged.bed')
    output:
        os.path.join(config['dir']['out'], 'results', 'bedtools.closest.out')
    shell:
        'bedtools closest -a {input.bed} -b {input.ref} -id -s -D ref> {output}'

rule find_matched:
    input:
        os.path.join(config['dir']['out'], 'results', 'bedtools.closest.out')
    output:
        os.path.join(config['dir']['out'], 'results',
                    'matched.bedtools.closest.out')
    shell:
        "awk '{{if ($11 !~ /\./ ) {{ print }} }}' {input} > {output}"

rule find_unmatched:
    input:
        os.path.join(config['dir']['out'], 'results', 'bedtools.closest.out')
    output:
        os.path.join(config['dir']['out'], 'results',
                    'unmatched.bedtools.closest.out')
    shell:
        "awk '{{if ($11 ~ /\./) {{ print }} }}' {input} > {output}"

rule filter_matched:
    input:
        os.path.join(config['dir']['out'], 'results',
                    'matched.bedtools.closest.out')
    output:
        os.path.join(config['dir']['out'], 'results',
                    'filtered.matched.bedtools.closest.out')
    params:
        kb=config['params']['kb']
    # look for distances > -kb because upstream exons in genome will be reported
    # as negative distances
    shell:
        "awk '{{if ($NF > -{params.kb}) {{ print }} }}' {input} > {output}"

