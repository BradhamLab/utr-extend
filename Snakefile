import re
from glob import glob
import os
import utils

configfile: 'config.yaml'
shell.prefix('source activate alignment; ')

lib_match = re.compile(config['params']['library_regex'])
bio_reads = config['params']['cdna_match']


FASTQS = [x for x in glob(os.path.join(config['dir']['data']\
          + f'*{bio_reads}*'))]

fastq_in = {'.fastq': 'cat',
           '.gz': 'zcat'}

read_fastq = fastq_in[os.path.splitext(FASTQS[0])[-1]]
library_dict = {}
LIBRARIES = []
for each in FASTQS:
    library = lib_match.search(os.path.basename(each)).group(0)
    library_dict[library] = each
    LIBRARIES.append(library)

msg = "\nIdentified libraries: {}\n\n".format(", ".join(LIBRARIES))\
    + "Expanding 3' UTR regions using reads from:\n\t{}.\n\n".format(
      "\n\t".join(FASTQS)) + "Fastq read command: {}\n".format(read_fastq) 
print(msg)

def fastq_input(wildcards):
    fastq = library_dict[wildcards.library]
    return "{library}".join(lib_match.split(fastq))

rule all:
    input:
        os.path.join(config['dir']['out'], 'results',
                     'updated_annotations.gff')

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
        overhang=config['params']['readlength'] - 1
    output:
        os.path.join(config['STAR']['index'], 'Genome')
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.overhang} --genomeChrBinNbits {params.chr_n_bits}"

rule align_3_prime_utr_reads:
    input:
        index=os.path.join(config['STAR']['index'], 'Genome'),
        fastq=lambda wildcards: fastq_input(wildcards)
    output:
        bam=os.path.join(config['dir']['out'], "STAR", "{library}",
                         "Aligned.sortedByCoord.out.bam")
    params:
        index=os.path.join(config['STAR']['index']),
        prefix=os.path.join(config['dir']['out'], "STAR", "{library}") + '/',
        read=read_fastq
    shell:
        'STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate '
        '--readFilesCommand {params.read} --genomeDir {params.index} --outFileNamePrefix '
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

# Do not ignore overlaps -- this will not 'filter out' reads already aligned to
# exons; instead it will report the upstream exon. Reads already mapped to exons
# should be removed during filtering
rule find_closest:
    input:
        ref=os.path.join(config['dir']['out'], 'annotations',
                         'exon.annotations.bed'),
        bed=os.path.join(os.path.join(config['dir']['out'], 'bed', '{library}',
                         'sorted.aligned.bed'))
    output:
        temp(os.path.join(config['dir']['out'], 'results', '{library}',
                          'bedtools.closest.out'))
    shell:
        "bedtools closest -a {input.bed} -b {input.ref} -D b -iu -s > {output}"

# remove alignments with zero distance -- already aligned to exon
# remove alignment more than maximum distance away
# remove alignments where | align.start - exon.start | > max dist and
# vice versa for end coordinates -- required for spliced alignments. 
rule filter_closest:
    input:
        os.path.join(config['dir']['out'], 'results', '{library}',
                        'bedtools.closest.out')
    output:
        temp(os.path.join(config['dir']['out'], 'results', '{library}',
                          'filtered.closest.out'))
    params:
        bp=config['params']['bp']
    shell:
        "awk '{{if ($10 !~ /\./ && $NF <= {params.bp} && $NF > 0 && "
        "sqrt(($2 - $8)^2) < {params.bp} && sqrt(($3 - $9)^2) < {params.bp} ) "
        "{{ print }} }}' {input} > {output}"

rule combine_libraries:
    input:
        expand(os.path.join(config['dir']['out'], 'results', '{library}',
                     'filtered.closest.out'), library=LIBRARIES)
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

rule update_gff:
    input:
        gff=config['genome']['gff'],
        extensions=os.path.join(config['dir']['out'], 'results',
                                '3_prime_alignments.out')
    output:
        gff=os.path.join(config['dir']['out'], 'results', 'updated_annotations.gff')
    params:
        bp=config['params']['bp']
    script:
        "scripts/annotate_3prime_UTR.py"


# clean gff