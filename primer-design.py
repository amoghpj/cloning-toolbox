import pandas as pd
from Bio import SeqIO
from optparse import OptionParser
# import primer3 as p3
# import matplotlib.pyplot as plt

def get_options():
    parse = OptionParser()
    parse.add_option("-i", "--input", help="CSV file with gene names")
    opts, vals = parse.parse_args()
    return(opts)

def get_sequence(sequences, gene):
    for fasta in sequences:
        sgdid, genename, *rest = fasta.description.split(" ")
        if gene.upper() == genename:
            fullorf = fasta.seq
    return(fullorf)

def main(opts):
    print(f"Input file: {opts.input}")
    genelist = pd.read_csv(opts.input, header=None)
    mapdf = pd.read_csv("./data/sgd-genename.csv",header=None, names=["systematic", "gene"])
    mapdf.set_index("gene",inplace=True)
    sequences = SeqIO.parse("./data/orf_genomic_1000.fasta",'fasta')
    ## Loop over each row of input
    for i, gene in enumerate(genelist[0].values):
        systematic = mapdf.loc[gene.upper()].systematic
        print(i, gene, systematic)
        full = get_sequence(sequences, gene)
        ## 1. Get donor DNA
        donor_f, donor_r = get_donor(full)
        ## 2. Get guides
        primer_list = get_primers(full)
        ## 3. Get PCR primers
        guides = get_guides(full)

def get_donor(full):
    """
    Not implemented
    """
    donor_f, donor_r = "", ""
    return(donor_r, donor_r)
def get_primers(full):
    """
    Not implemented
    """
    return([])

def get_guides(full):
    """
    Not implemented
    """
    return([])

if __name__ == '__main__':
    opts = get_options()
    main(opts)

# seq = # 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAA'
#   # TGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
# o = p3.designPrimers(
#     {
#         'SEQUENCE_ID': 'MH1000',
#         'SEQUENCE_TEMPLATE': seq,
#         'SEQUENCE_INCLUDED_REGION': [0, len(seq)]
#     },
#     {
#         'PRIMER_OPT_SIZE': 20,
#         'PRIMER_PICK_INTERNAL_OLIGO': 1,
#         'PRIMER_INTERNAL_MAX_SELF_END': 8,
#         'PRIMER_MIN_SIZE': 18,
#         'PRIMER_MAX_SIZE': 25,
#         'PRIMER_OPT_TM': 60.0,
#         'PRIMER_MIN_TM': 57.0,
#         'PRIMER_MAX_TM': 63.0,
#         'PRIMER_MIN_GC': 20.0,
#         'PRIMER_MAX_GC': 80.0,
#         'PRIMER_MAX_POLY_X': 100,
#         'PRIMER_INTERNAL_MAX_POLY_X': 100,
#         'PRIMER_SALT_MONOVALENT': 50.0,
#         'PRIMER_DNA_CONC': 50.0,
#         'PRIMER_MAX_NS_ACCEPTED': 0,
#         'PRIMER_MAX_SELF_ANY': 12,
#         'PRIMER_MAX_SELF_END': 8,
#         'PRIMER_PAIR_MAX_COMPL_ANY': 12,
#         'PRIMER_PAIR_MAX_COMPL_END': 8,
#         'PRIMER_PRODUCT_SIZE_RANGE': [[100,min(len(seq),900)]],
#     })

# from dna_features_viewer import GraphicFeature, GraphicRecord
# features=[GraphicFeature(start=0, end=len(seq), strand=+1, color="#000000",
#                    label="Template")
# ]
# for i in range(0,5):
#     start, end = o[f"PRIMER_LEFT_{i}"]
#     features.append(GraphicFeature(start=start, end=start+20, strand=1, color="#ffd700",
#                    label=o[f'PRIMER_LEFT_{i}_SEQUENCE']))
#     start, end = o[f"PRIMER_RIGHT_{i}"]
#     features.append(GraphicFeature(start=start, end=start+20, strand=-1, color="#ffd700",
#                    label=o[f'PRIMER_RIGHT_{i}_SEQUENCE']))

# record = GraphicRecord(sequence_length=len(seq)+100 , features=features)
# record.plot(figure_width=10)
# plt.show()


