"""
Description: Script to generate sgRNA candidates, donor sequences, and colony verification primers.
Author: Amogh Jalihal
Date: 2022-05-08
"""
import os
import subprocess
import primer3 as p3
import pandas as pd
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import BiopythonTranslator
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord

def get_options():
    parse = OptionParser()
    parse.add_option("-i", "--input", help="CSV file with gene names")
    opts, vals = parse.parse_args()
    return(opts)

def get_donor(full, gene):
    """
    Get donor sequencs
    """
    arm_size = 50
    overlap_size = 20
    donor_f, donor_r = "", ""
    promoter_arm = full[1000-arm_size:1000]
    terminator_arm = full[-1000:-1000 + arm_size].reverse_complement()
    overlap_l = promoter_arm[int(-overlap_size/2):]
    overlap_r = terminator_arm[int(-overlap_size/2):]
    donor_f = promoter_arm + overlap_r.reverse_complement()
    donor_r = terminator_arm + overlap_l.reverse_complement()
    seqdict = {}

    features=[SeqFeature(FeatureLocation(1000-arm_size, 1000), 
                         type = "promoter"),
              SeqFeature(FeatureLocation(1000, len(full)-1000), 
                         type = "CDS"),
              SeqFeature(FeatureLocation(len(full)-1000, 
                                         len(full) - 1000 + arm_size),
                         type = "terminator")]
    seqdict[f"{gene}-donor-f"] = donor_f
    seqdict[f"{gene}-donor-r"] = donor_r
    return(seqdict,features)

def get_primers(full, gene):
    """
    Return top PCR primers.
    """
    p3_global_parameters = {
        'PRIMER_OPT_SIZE': 20,
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
    CDS_END = CDS_START + min(CDS_LENGTH, 1000)
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
    seqdict[f"{gene}-col-f"] = f_primer_seq
    wt_primers = p3.designPrimers({
        'SEQUENCE_ID': f"{gene}_wt",
        'SEQUENCE_TEMPLATE': str(wt),
        "PRIMER_LEFT": f_primer_seq,
        "PRIMER_TARGET": str(wt)[-100:],
        'PRIMER_PRODUCT_SIZE_RANGE': [[800, len(wt)]]},
                                  p3_global_parameters)
    del_primers = p3.designPrimers({
        'SEQUENCE_ID': f"{gene}_del",
        'SEQUENCE_TEMPLATE': str(deletion),
        "PRIMER_LEFT": f_primer_seq,
        'PRIMER_PRODUCT_SIZE_RANGE': [[100,400]]}, 
                         p3_global_parameters)
    pr_loc = full.find(f_primer_seq)
    for i in range(4):
        start, end = wt_primers[f"PRIMER_RIGHT_{i}"]
        loc = full.find(Seq(wt_primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]).reverse_complement())
        if (loc - pr_loc > 600) and (loc - pr_loc < 1000):
            break
        
    rec.append(SeqFeature(FeatureLocation(PROMOTER_START + start -end, 
                                          PROMOTER_START + start ),  
                          type = "pwt_r", strand=-1))
    wt_loc = full.find(Seq(wt_primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]).reverse_complement())
    seqdict[f"{gene}-col-r1({wt_loc - pr_loc})"] = wt_primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]

    i = 0
    start, end = del_primers[f"PRIMER_RIGHT_{i}"]
    rec.append(SeqFeature(FeatureLocation(TERM_END - start, TERM_END - start + end), 
                          type = "pdel_r", strand=-1))
    del_loc = full.find(Seq(del_primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]).reverse_complement())
    print(pr_loc)
    print(wt_loc)
    print(del_loc)

    seqdict[f"{gene}-col-r2({del_loc - CDS_LENGTH - pr_loc })"] = del_primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    return(seqdict, rec)

def get_guides(full,gene):
    """
    Call CCTop to generate candidate guide sequences
    """
    seqdict = {} 
    rec = []
    INDEXPATH = os.path.expanduser("~") + "/s_cerevisiae/s_cerevisiae"
    if not os.path.exists(f"{gene}-dir"):
        os.mkdir(f"{gene}-dir")

    with open(f"{gene}-dir/{gene}-sequence.fasta", "w") as outfile:
        outfile.write(f">{gene}\n")
        outfile.write(f"{full[1000:-1000]}")
    p = subprocess.Popen(["cctop","--input", f"{gene}-dir/{gene}-sequence.fasta",
                      "--index", INDEXPATH,
                       "--output", f"{gene}-dir"])
    p.wait()
    f = [_f for _f in glob.glob(f"{gene}-dir/*.fasta") if "-sequence" not in _f][0]
    with open(f, "r") as infile:
        for i in range(6):
            header, sequence = infile.readline(), infile.readline()
            name = header.strip().replace(">","")
            gibson_construct = Seq("GACTTT" + sequence.strip()[:-3] + "GTTT")
            seqdict[f"{gene}-sgrna{name}-f"] = gibson_construct[:-4]
            seqdict[f"{gene}-sgrna{name}-r"] = gibson_construct[4:].reverse_complement()
            start = full.find(sequence.strip())
            rec.append(SeqFeature(FeatureLocation(start, start + len(sequence)), 
                         type = f"guide-{i}"))
    return(seqdict, rec)

def get_sequence(gene):
    sequences = SeqIO.parse("./data/orf_genomic_1000.fasta",'fasta')
    for fasta in sequences:
        sgdid, genename, *rest = fasta.description.split(" ")
        if gene.upper() == genename:
            fullorf = fasta.seq
            break
    rec = SeqRecord(fullorf, id=gene, annotations={"molecule_type":"DNA",
                                                "organism":"S. cerevisiae"},
                                                features=[])
    return(fullorf, rec)

def main(opts):
    print(f"Input file: {opts.input}")
    genelist = pd.read_csv(opts.input, header=None)
    mapdf = pd.read_csv("./data/sgd-genename.csv",header=None, names=["systematic", "gene"])
    mapdf.set_index("gene",inplace=True)

    ## Loop over each row of input
    for i, gene in enumerate(genelist[0].values):
        alloligos = {}
        systematic = mapdf.loc[gene.upper()].systematic
        print(i, gene, systematic)
        full, record = get_sequence(gene)
        
        ## 1. Get donor DNA
        print("Designing donors...")
        seqdict, rec = get_donor(full, gene)
        alloligos.update(seqdict)
        record.features.extend(rec)

        # ## 2. Get guides
        # print("Designing guides...")
        # seqdict, rec = get_guides(full, gene)
        # alloligos.update(seqdict)
        # record.features.extend(rec)

        ## 3. Get PCR primers
        print("Designing colony pcr primers...")
        seqdict, rec = get_primers(full, gene)
        alloligos.update(seqdict)
        record.features.extend(rec)

        print("...done!")

        print("Writing genbank file...")
        SeqIO.write(record, f"{gene}.gb", "genbank")
        print("Writing fasta sequence list...")
        with open(f"{gene}.fasta", "w") as outfile:
            for k, val in alloligos.items():
                outfile.write(f">{k}\n")
                outfile.write(f"{val}\n")

        print("Writing to csv...")
        with open(f"{gene}.csv", "w") as outfile:
            for k, val in alloligos.items():
                outfile.write(f"{k},{val}\n")

        print("Generating graphic...")
        graphic_record = BiopythonTranslator().translate_record(f"{gene}.gb")
        graphic_record = graphic_record.crop((400,len(full)-400 ))
        graphic_record.plot(figure_width=6, strand_in_label_threshold=7)
        plt.savefig(f"{gene}.png")

if __name__ == '__main__':
    opts = get_options()
    main(opts)



