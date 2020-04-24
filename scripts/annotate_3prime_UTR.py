from collections import namedtuple
import warnings

import gffutils
import re



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


def update_lengths(db, extensions, max_dist=5000):
    utr_to_mrna = dict()
    update_features = []
    # new_features = []
    delete_features = []
    with open(extensions, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            align = Alignment(chr=line_split[0],
                            start=int(line_split[1]) + 1, # bed to gff coords
                            end=int(line_split[2]) + 1,
                            strand=line_split[4],
                            exon=line_split[-1])
            exon = db[align.exon]
            # get gene
            gene = next(db.parents(exon, featuretype='gene'))
            mrna = next(db.parents(exon, featuretype='mRNA'))
            # if exon has already been extended to alignment length, no update nec
            if exon.end >= align.end and exon.start <= align.start:
                print('already expanded')
                continue

            utr = gffutils.Feature(seqid=exon.chrom,
                                   source='utr-extend',
                                   featuretype='three_prime_UTR',
                                   start=exon.end,
                                   end=align.end,
                                   id="{}:{}".format(exon.id, 'three_prime_UTR'),
                                   strand=exon.strand)
            
            # ensure extensions don't exceed a maximum distance
            if align.strand == '+' and align.end - align.start > max_dist:
                utr.end = exon.end + max_dist 

            elif align.strand == '-':
                if align.end - align.start > max_dist:
                    align.start = exon.start - max_dist
                utr.start = align.start
                utr.end = exon.start
            # check for overlapping utr annotations in the given gene
            present_utrs = [x for x in db.region(exon,
                                                 featuretype='three_prime_UTR',
                                                 strand=exon.start)\
                            if next(db.parents(x, featuretype='mRNA')).id == mrna.id]
            # if current utr in exon, set new start/stop. Set old id to be deleted
            if len(present_utrs) > 0:
                sources = set()
                for each in present_utrs:
                    if each.start < utr.start:
                        utr.start = each.start
                        sources.add('utr-extend')
                    if each.stop > utr.stop:
                        utr.stop = each.stop
                        sources.add('utr-extend')
                    # fine to delete even if we take id since it is added back
                    # in with proceeding update() call -- probably
                    delete_features.append(each)
                    sources.add(each.source)
                if len(present_utrs) == 1:
                    utr.id = present_utrs[0].id
                utr.source = ','.join(sorted(sources))

            # increase related boundaries
            if utr.strand == '+':
                if utr.end > exon.end:
                    exon.end = utr.end
                    if 'utr-extend' not in exon.source:
                        exon.source += ',utr-extend'
                if utr.end > gene.end:
                    gene.end = utr.end
                if utr.end > mrna.end:
                    mrna.end = utr.end
            elif utr.strand == '-':
                if utr.start < exon.start:
                    exon.start = utr.start
                    if 'utr-extend' not in exon.source:
                        exon.source += ',utr-extend'
                if utr.start < gene.start:
                    gene.start = utr.start
                if utr.start < mrna.start:
                    mrna.start = utr.start

            # check overlapping exons with same parent gene
            overlapping_exons = [x for x in db.region(exon, featuretype='exon',
                                                      strand=exon.strand)\
                                 if x.id != exon.id and\
                                 next(db.parents(x, featuretype='gene')).id == gene.id]
            if len(overlapping_exons) > 0:
                warnings.warn("Overlapping exons! Current expanded UTR:\n\t"\
                              "{}.\n\nOverlapping Exons:\n\t{}".format(
                               utr, "\n\t".join([str(x) for x in overlapping_exons])))
            utr.attributes['ID'] = utr.id
            utr.attributes['Parent'] = mrna.id
            # new_features.append(utr)
            update_features += [gene, exon, mrna, utr]
            utr_to_mrna[utr] = mrna
            

    db = db.delete(delete_features)
    # db = db.update(new_features, merge_strategy='create_unique')
    db = db.update(update_features, merge_strategy='replace')

    for utr, mrna in utr_to_mrna.items():
        # try to set relation from gene to utr -- fails if relationship is
        # present so just pass
        try:
            db = db.add_relation(mrna, utr, level=1)
        except:
            pass
    return db



if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        fn = "/projectnb/bradham/data/ReferenceSequences/wray-genome/L_var_clean.gff"
        db = gffutils.create_db(snakemake.input['gff'], ':memory:', 
                                merge_strategy='create_unique')
        db = update_lengths(db, snakemake.input['extensions'],
                            snakemake.params['bp'])
        with open(snakemake.output['gff'], 'w') as fout:
            for gene in db.features_of_type('gene', order_by=('seqid', 'start')):
                fout.write(str(gene) + '\n')
                for f in db.children(gene, order_by='start'):
                    fout.write(str(f) + '\n')