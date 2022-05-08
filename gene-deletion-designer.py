import os
import subprocess
import primer3 as p3
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import BiopythonTranslator
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord
# features=[GraphicFeature(start=0, end=len(seq), strand=+1, color="#000000",
#                   label="Template")]
# 


def get_options():
    parse = OptionParser()
    parse.add_option("-i", "--input", help="CSV file with gene names")
    opts, vals = parse.parse_args()
    return(opts)

def get_donor(full, gene):
    """
    Not implemented
    """
    arm_size = 50
    overlap_size = 20
    donor_f, donor_r = "", ""
    promoter_arm = full[1000-arm_size:1000]
    terminator_arm = full[-1000:-1000 + arm_size].reverse_complement()
    overlap = promoter_arm[-overlap_size]
    donor_f = promoter_arm + overlap
    donor_r = terminator_arm + Seq(overlap).reverse_complement()
    seqdict = {}
    rec = SeqRecord(full, id=gene, annotations={"molecule_type":"DNA", 
                                                "organism":"S. cerevisiae",
                                                "comments": "Generated automatically"},
                    features=[SeqFeature(FeatureLocation(1000-arm_size, 1000), 
                                         type = "promoter"),
                              SeqFeature(FeatureLocation(1000, len(full)-1000), 
                                         type = "CDS"),
                              SeqFeature(FeatureLocation(len(full)-1000, 
                                                         len(full) - 1000 + arm_size),
                                         type = "terminator")])
    seqdict[f"{gene}_donor_f"] = donor_f
    seqdict[f"{gene}_donor_r"] = donor_r
    return(seqdict, rec)

def get_primers(full, gene):
    """
    Return top PCR primers.
    """
    p3_global_parameters = {
        'PRIMER_OPT_SIZE': 20,
        # 'PRIMER_PICK_INTERNAL_OLIGO': 1,
        # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
        # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    }
    PROMOTER_START = 700
    CDS_LENGTH = len(full) - 2000
    CDS_START = 1000
    CDS_END = CDS_START + CDS_LENGTH
    TERM_END = 1000 + CDS_LENGTH + 500
    TERM_START = len(full) - 1000 

    seqdict = {}
    rec = []

    ### Define semantic sequence stretches
    wt = full[PROMOTER_START: CDS_END]
    deletion = full[TERM_START :TERM_END]
    promoter = full[PROMOTER_START:CDS_START]

    ### Get forward primers in promoter region.
    ## These will be used to design reverse primers further down
    force_forward =  p3.designPrimers({
        'SEQUENCE_ID': f"{gene}_pr",
        'SEQUENCE_TEMPLATE': str(promoter),
        'PRIMER_PRODUCT_SIZE_RANGE': [[100,200]]}, 
                         p3_global_parameters)

    f_primer_seq = force_forward["PRIMER_LEFT_1_SEQUENCE"]
    start, end = force_forward[f"PRIMER_LEFT_1"]
    rec.append(SeqFeature(FeatureLocation(PROMOTER_START + start, PROMOTER_START  + start + end), 
                                     type = "pwt_f"))
    seqdict[f"{gene}_f"] = f_primer_seq
    wt_primers = p3.designPrimers({
        'SEQUENCE_ID': f"{gene}_wt",
        'SEQUENCE_TEMPLATE': str(wt),
        "PRIMER_LEFT": f_primer_seq,
        "PRIMER_TARGET": str(wt)[:-100],
        'PRIMER_PRODUCT_SIZE_RANGE': [[700, 900]]},
        #'SEQUENCE_INCLUDED_REGION': [0, len(wt)]}, 
                                  p3_global_parameters)
    del_primers = p3.designPrimers({
        'SEQUENCE_ID': f"{gene}_del",
        'SEQUENCE_TEMPLATE': str(deletion),
        "PRIMER_LEFT": f_primer_seq,
        'PRIMER_PRODUCT_SIZE_RANGE': [[100,400]]}, 
                         p3_global_parameters)
    
    i = 1
    start, end = wt_primers[f"PRIMER_RIGHT_{i}"]
    rec.append(SeqFeature(FeatureLocation(PROMOTER_START + start -end, 
                                          PROMOTER_START + start ),  
                                    type = "pwt_r", strand=-1))
    seqdict[f"{gene}_wt_r"] = wt_primers[f"PRIMER_RIGHT_{i}"]

    start, end = del_primers[f"PRIMER_RIGHT_{i}"]
    rec.append(SeqFeature(FeatureLocation(TERM_END - start, TERM_END - start + end), 
                          type = "pdel_r", strand=-1))
    seqdict[f"{gene}_del_r"] = del_primers[f"PRIMER_RIGHT_{i}"]
    return(seqdict, rec)

def get_guides(full,gene):
    """
    Call CCTop to generate candidate guide sequences
    """
    seqdict = {} 
    rec = []
    INDEXPATH = os.path.expanduser("~") + "/s_cerevisiae/s_cerevisiae"
    os.mkdir(f"{gene}-dir")
    with open(f"{gene}-dir/{gene}-sequence.fasta", "w") as outfile:
        outfile.write(f"> {gene}\n")
        outfile.write(f"{full}")
    p = subprocess.Popen(["cctop","--input", f"{gene}-dir/{gene}-sequence.fasta",
                      "--index", INDEXPATH,
                       "--output", f"{gene}-dir"])
    p.wait()
    f = glob.glob(f"{gene}-dir/*.fasta")[0]
    with open(f, "r") as infile:
        for line in infile.readlines()[:5]:
            print(line)
    return(seqdict, rec)

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
    alloligos = {}
    ## Loop over each row of input
    for i, gene in enumerate(genelist[0].values):
        systematic = mapdf.loc[gene.upper()].systematic
        print(i, gene, systematic)
        full = get_sequence(sequences, gene)
        
        ## 1. Get donor DNA
        seqdict, rec = get_donor(full, gene)
        alloligos.update(seqdict)

        ## 2. Get guides
        seqdict, rec = get_primers(full, gene)
        alloligos.update(seqdict)

        ## 3. Get PCR primers
        seqdict, rec = get_guides(full)
        alloligos.update(seqdict)

        rec.features.extend(primer_list)
        SeqIO.write(rec, f"{gene}.gb", "genbank")
        with open(f"{gene}.fasta", "w") as outfile:
            for k, val in alloligos.items():
                outfile.write(f"{k}\n")
                outfile.write(f"{val}\n")
        graphic_record = BiopythonTranslator().translate_record(f"{gene}.gb")
        fig = plt.figure(figsize=(10,5))
        axes = [fig.add_subplot(2,1,i) for i in range(1,3)]
        graphic_record = graphic_record.crop((400,len(full) ))
        graphic_record.plot(strand_in_label_threshold=7, ax=axes[0])
        cropped = graphic_record.crop((950,1050))
        cropped.plot(ax=axes[1])
        cropped.plot_sequence(ax=axes[1])
        plt.savefig(f"{gene}.png")

if __name__ == '__main__':
    opts = get_options()
    main(opts)



