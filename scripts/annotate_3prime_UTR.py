from collections import namedtuple

import gffutils
import re

# TODO account for gff/bed coordinate offsets

fn = "/projectnb/bradham/data/ReferenceSequences/wray-genome/L_var_clean.gff"
db = gffutils.create_db(fn, ':memory:', merge_strategy='create_unique')

# get transcript
# get exon
# extend exon end to alignment end
# create three_prime_UTR from exon end to 
Alignment = namedtuple('Alignment',
                       ['chr',
                        'start',
                        'end',
                        'strand',
                        'exon'])


def update_lengths(db, line, max_dist=5000):
    # this is assuming '+' strand 
    line_split = line.strip().split('\t')
    align = Alignment(chr=line_split[0],
                      start=int(line_split[1]) + 1, # bed to gff coords
                      end=int(line_split[2]),
                      strand=line_split[4],
                      exon=line_split[-1])
    exon = db[align.exon]
    # get gene
    gene = next(db.parents(exon.id, featuretype='gene'))
    exon_name = re.search(gene.id + '[^:]*', exon.id).group(0)
    # if exon has already been extended to alignment length, no update nec
    if exon.end >= align.end and exon.start <= align.start:
        print('already expanded')
        return db

    if align.strand == '+':
        # ensure extensions don't exceed a maximum distance
        if align.end - align.start > max_dist:
            align.end = exon.end + max_dist 
        utr = gffutils.Feature(seqid=f"{exon_name}:three_prime_utr",
                               source='utr-extend',
                               featuretype='three_prime_UTR',
                               start=exon.end,
                               end=align.end,
                               strand=exon.strand)
    elif align.strand == '-':
        if align.end - align.start > max_dist:
            align.start = exon.start - max_dist
        utr = gffutils.Feature(seqid=f"{exon_name}:three_prime_utr",
                               source='utr-extend',
                               featuretype='three_prime_UTR',
                               start=align.start,
                               end=exon.start,
                               strand=exon.strand)
 
    # iterate through current utrs if overlap in exon consider the same
    present_utrs = db.children(gene.id, featuretype='three_prime_UTR')
    n_utrs = 0
    new_utr = True
    for each in present_utrs:
        n_utrs += 1
        # check if it's contained in current exon
        if exon.strand == '+':
            if exon.start <= each.start and each.start <= exon.end:
                if utr.end > each.end:
                    each.end = utr.end
                if utr.start < each.start:
                    each.start = utr.start
                utr = each
                new_utr = False
                break
        elif exon.strand == '-':
            if exon.start <= each.end and each.end <= exon.end:
                if utr.end > each.end:
                    each.end = utr.end
                if utr.start < each.start:
                    each.start = utr.start
                new_utr = False
                break


    if new_utr and n_utrs > 0:
        utr.id += f"_{n_utrs - 1}"
    # + strand
    if utr.end > gene.end:
        gene.end = utr.end
    if utr.end > exon.end:
        exon.end = utr.end
    # - strand stuff
    if utr.start < gene.start:
        gene.start = utr.start
    if utr.start < exon.start:
        exon.start = utr.start

    if new_utr:
        db = db.add_relation(gene, utr, level=3)

    db = db.update([gene, exon, utr], merge_strategy='replace')
    
    return db