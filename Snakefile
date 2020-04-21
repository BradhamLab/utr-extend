import re
from glob import glob
import os
import utils

configfile: 'config.yaml'
shell.prefix('source activate alignment; ')

lib_match = re.compile(config['params']['library_regex'])
bio_reads = config['params']['cdna_match']
# change this
FASTQS = [x for x in glob(os.path.join(config['dir']['data']\
          + f'*{bio_reads}*'))]
library_dict = {}
LIBRARIES = []
for each in FASTQS:
    library = lib_match.search(os.path.basename(each)).group(0)
    library_dict[library] = each
    LIBRARIES.append(library)
print(LIBRARIES)

def fastq_input(wildcards):
    fastq = library_dict[wildcards.library]
    return "{library}".join(lib_match.split(fastq))

rule all:
    input:
        os.path.join(config['dir']['out'], 'results',
                     '3_prime_alignments.out')

rule gff_to_bed:
    input:
        gtf=config['genome']['gff']
    output:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'annotations.bed')
    shell:
        'gff2bed < {input.gtf} > {output.bed}'

rule sort_annotations:
    input:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'annotations.bed')
    output:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'sorted.annotations.bed')
    shell:
        'sort -k1,1 -k2,2n {input.bed} > {output.bed}'

rule extract_exon_annotations:
    input:
        bed=os.path.join(config['dir']['out'], 'annotations',
                         'sorted.annotations.bed')
    output:
        bed=os.path.join(config['dir']['out'], 'annotations',
                'exon.annotations.bed')
    shell:
        "awk '{{if ($8 ~ /exon/) {{ print }} }}' {input} > {output}"

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

rule align_3_prime_utr_reads:
    input:
        index=os.path.join(config['STAR']['index'], 'Genome'),
        fastq=lambda wildcards: fastq_input(wildcards)
    output:
        bam=os.path.join(config['dir']['out'], "STAR", "{library}",
                         "Aligned.sortedByCoord.out.bam")
    params:
        index=os.path.join(config['STAR']['index']),
        prefix=os.path.join(config['dir']['out'], "STAR", "{library}") + '/'
    shell:
        'STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate '
        '--readFilesCommand zcat --genomeDir {params.index} --outFileNamePrefix '
        '{params.prefix} --readFilesIn {input.fastq} --outReadsUnmapped FASTX'

# switch sam to sorted bam
rule bam_to_bed:
    input:
        bam=os.path.join(config['dir']['out'], "STAR",
                         "{library}", "Aligned.sortedByCoord.out.bam")
    output:
        bed=temp(os.path.join(config['dir']['out'], 'bed', '{library}',
                              'aligned.bed'))
    shell:
        'bedtools bamtobed -i {input.bam} > {output.bed}'

rule sort_alignment_bed:
    input:
        bed=os.path.join(config['dir']['out'], 'bed', '{library}',
                         'aligned.bed')
    output:
        bed=temp(os.path.join(config['dir']['out'], 'bed', '{library}',
                         'sorted.aligned.bed'))
    shell:
        'sort -k1,1 -k2,2n {input.bed} > {output.bed}'
    
rule find_closest:
    input:
        ref=os.path.join(config['dir']['out'], 'annotations',
                         'exon.annotations.bed'),
        bed=os.path.join(os.path.join(config['dir']['out'], 'bed', '{library}',
                         'sorted.aligned.bed'))
    output:
        temp(os.path.join(config['dir']['out'], 'results', '{library}',
                    'bedtools.closet.out'))
    shell:
        "bedtools closest -a {input.bed} -b {input.ref} -D b -io -iu -s > {output}"

rule filter_closest:
    input:
        os.path.join(config['dir']['out'], 'results', '{library}',
                        'bedtools.closet.out')
    output:
        temp(os.path.join(config['dir']['out'], 'results', '{library}',
                     'merged.closest.out'))
    params:
        bp=config['params']['bp']
    shell:
        "awk '{{if ($10 !~ /\./ && $NF <= {params.bp}) {{ print }} }}' {input} > {output}"

rule combine_libraries:
    input:
        expand(os.path.join(config['dir']['out'], 'results', '{library}',
                     'merged.closest.out'), library=LIBRARIES)
    output:
        temp(os.path.join(config['dir']['out'], 'results', 'combined.out'))
    shell:
        "cat {input} > {output}"

# then merge
rule group_library_by_exon:
    input:
        os.path.join(config['dir']['out'], 'results', 'combined.out')
    output:
        os.path.join(config['dir']['out'], 'results', '3_prime_alignments.out')
    shell:
        "bedtools groupby -i {input} -g 10 -c 1,2,3,4,6,7,8,9,10 "
        "-o distinct,min,max,count_distinct,distinct,distinct,distinct,distinct,distinct "
        " | cut -f2- > {output}"