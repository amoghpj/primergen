import sys
from Bio import SeqIO, Seq
import os
import primer3 as p3
import pandas as pd

NUM_PRIMERS_PER_CDS = 5
NUM_PRIMERS_TO_GENERATE = 300

def get_primers(record):
    """
    Return top PCR primers.
    """
    p3_global_parameters = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_MIN_TM': 59.0,
        'PRIMER_MAX_TM': 64.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    }
    seqcollect = []
    rec = []

    primers =  p3.design_primers({
        'SEQUENCE_ID': f"{record.id}",
        'SEQUENCE_TEMPLATE': record.seq,
        'PRIMER_PRODUCT_SIZE_RANGE': [[90,120]]}, 
                         p3_global_parameters)
    ok = True
    i = 1
    while ok:
        try:
            # Forward
            startf, endf = primers[f"PRIMER_LEFT_{i}"]
            # Reverse
            startr, endr = primers[f"PRIMER_RIGHT_{i}"]
            seqcollect.append({"gene":record.id,
                               "primer_id":i,
                               "forward_start":startf,
                               "forward_end":startf + endf,
                               "reverse_start":startr,
                               "reverse_end":startr + endr,
                               "primer_pair_penalty":primers[f"PRIMER_PAIR_{i}_PENALTY"],
                               "forward_tm":primers[f"PRIMER_LEFT_{i}_TM"],
                               "reverse_tm":primers[f"PRIMER_RIGHT_{i}_TM"],
                               "forward_gc":primers[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                               "reverse_gc":primers[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                               "reverse_end":startr + endr,
                               "amplicon_sequence":record.seq[startf:endr]
                               "amplicon_size":abs(startr + endr - (startf)),
                               "product_tm":primers[f"PRIMER_PAIR_{i}_PRODUCT_TM"],
                               "forward_primer":primers[f"PRIMER_LEFT_{i}_SEQUENCE"],
                               "reverse_primer":primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]})
        except:
            ok = False
        i += 1
        if i > NUM_PRIMERS_PER_CDS:
            ok = False

    return(seqcollect)

def is_valid(s):
    isvalid = True
    if "ribosomal protein S9/S16" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein S8" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein L10" in s:
        isvalid = False
        print(f"removing {s}")
    if "grpe" in s:
        print(f"removing {s}")
        isvalid = False
    if "rimp N-terminal domain" in s:
        print(f"removing {s}")
        isvalid = False
    if "polyribonucleotide nucleotidyltransferase, RNA binding domain" in s:
        print(f"removing {s}")
        isvalid = False
    if "16S rRNA (cytosine(1402)-N(4))-methyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "peptide chain release factor 1" in s:
        print(f"removing {s}")
        isvalid = False
    if "peptide chain release factor 2" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bS20" in s:
        print(f"removing {s}")
        isvalid = False
    if "rRNA maturation RNase YbeY" in s:
        print(f"removing {s}")
        isvalid = False
    if "RIP metalloprotease RseP" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bL17" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bL21" in s:
        print(f"removing {s}")
        isvalid = False
    if "signal recognition particle-docking protein FtsY" in s:
        print(f"removing {s}")
        isvalid = False
    if "cell division protein FtsZ" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosome-binding factor A" in s:
        print(f"removing {s}")
        isvalid = False
    if "riboflavin biosynthesis protein RibF" in s:
        print(f"removing {s}")
        isvalid = False
    if "Holliday junction DNA helicase RuvA" in s:
        print(f"removing {s}")
        isvalid = False
    if "SsrA-binding protein" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA (guanine(37)-N(1))-methyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosome silencing factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "GTP-binding protein YchF" in s:
        print(f"removing {s}")
        isvalid = False
    if "16S rRNA (guanine(966)-N(2))-methyltransferase RsmD" in s:
        print(f"removing {s}")
        isvalid = False
    if "trigger factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "translation elongation factor Ts" in s:
        print(f"removing {s}")
        isvalid = False
    if "16S rRNA (guanine(527)-N(7))-methyltransferase RsmG" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bL9" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bS6" in s:
        print(f"removing {s}")
        isvalid = False
    if "translation initiation factor IF-3" in s:
        print(f"removing {s}")
        isvalid = False
    if "RNA methyltransferase, TrmH family, group 3" in s:
        print(f"removing {s}")
        isvalid = False
    if "excinuclease ABC subunit C" in s:
        print(f"removing {s}")
        isvalid = False
    if "putative transcription antitermination factor YqgF" in s:
        print(f"removing {s}")
        isvalid = False
    if "CTP synthase" in s:
        print(f"removing {s}")
        isvalid = False
    if "alanine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "chromosomal replication initiator protein DnaA" in s:
        print(f"removing {s}")
        isvalid = False
    if "ATP-dependent Clp protease, ATP-binding subunit ClpX" in s:
        print(f"removing {s}")
        isvalid = False
    if "isoleucine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "leucine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "methionine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "serine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA repair protein RadA" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA (5-methylaminomethyl-2-thiouridylate)-methyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA pseudouridine(55) synthase" in s:
        print(f"removing {s}")
        isvalid = False
    if "cysteine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "GTP-binding protein Era" in s:
        print(f"removing {s}")
        isvalid = False
    if "histidine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "phospho-N-acetylmuramoyl-pentapeptide-transferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "arginine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "aspartate--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "methionyl-tRNA formyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "phenylalanine--tRNA ligase, alpha subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "phenylalanine--tRNA ligase, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "translation initiation factor IF-2" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosome recycling factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "putative oxygen-independent coproporphyrinogen III oxidase" in s:
        print(f"removing {s}")
        isvalid = False
    if "transcription-repair coupling factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA polymerase I" in s:
        print(f"removing {s}")
        isvalid = False
    if "recombination protein RecR" in s:
        print(f"removing {s}")
        isvalid = False
    if "eifnuclease ABC subunit B" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA repif protein RecN" in s:
        print(f"removing {s}")
        isvalid = False
    if "Holliday junction DNA helicase RuvB" in s:
        print(f"removing {s}")
        isvalid = False
    if "ATP-dependent DNA helicase RecG" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA polymerase III, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bS1" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal RNA small subunit methyltransferase A" in s:
        print(f"removing {s}")
        isvalid = False
    if "preprotein translocase, SecG subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "transcription termination/antitermination factor NusG" in s:
        print(f"removing {s}")
        isvalid = False
    if "adenylosuccinate lyase" in s:
        print(f"removing {s}")
        isvalid = False
    if "signal recognition particle protein" in s:
        print(f"removing {s}")
        isvalid = False
    if "preprotein translocase, SecA subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "preprotein translocase, SecE subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "preprotein translocase, SecY subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS3" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS2" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS4" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS5" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS7" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein bL20" in s:
        print(f"removing {s}")
        isvalid = False
    if "ATP synthase F1, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL22" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA gyrase, B subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA gyrase, A subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL13" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL15" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL24" in s:
        print(f"removing {s}")
        isvalid = False
    if "UDP-N-acetylmuramate--L-alanine ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "UDP-N-acetylmuramoylalanine--D-glutamate ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA polymerase III, delta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ATP synthase F1, gamma subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL16" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL1" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL2" in s:
        print(f"removing {s}")
        isvalid = False
    if "inosine-5'-monophosphate dehydrogenase" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA primase" in s:
        print(f"removing {s}")
        isvalid = False
    if "elongation factor 4" in s:
        print(f"removing {s}")
        isvalid = False
    if "GTP-binding protein TypA/BipA" in s:
        print(f"removing {s}")
        isvalid = False
    if "pantetheine-phosphate adenylyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL11" in s:
        print(f"removing {s}")
        isvalid = False
    if "transcription antitermination factor NusB" in s:
        print(f"removing {s}")
        isvalid = False
    if "transcription termination factor NusA" in s:
        print(f"removing {s}")
        isvalid = False
    if "protein RecA" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA-directed RNA polymerase, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA-directed RNA polymerase, alpha subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "UMP kinase" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribonuclease III" in s:
        print(f"removing {s}")
        isvalid = False
    if "16S rRNA processing protein RimM" in s:
        print(f"removing {s}")
        isvalid = False
    if "chaperone protein DnaK" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA-directed RNA polymerase, beta' subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "DNA polymerase III, subunit gamma and tau" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA(Ile)-lysidine synthetase" in s:
        print(f"removing {s}")
        isvalid = False
    if "Obg family GTPase CgtA" in s:
        print(f"removing {s}")
        isvalid = False
    if "guanylate kinase" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosome-associated GTPase EngA" in s:
        print(f"removing {s}")
        isvalid = False
    if "50S ribosomal protein uL3" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uS11" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribosomal protein uL6" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA threonylcarbamoyl adenosine modification protein TsaD" in s:
        print(f"removing {s}")
        isvalid = False
    if "tRNA threonylcarbamoyl adenosine modification protein YeaZ" in s:
        print(f"removing {s}")
        isvalid = False
    if "50S ribosomal protein uL4" in s:
        print(f"removing {s}")
        isvalid = False
    if "hypothetical protein" in s:
        isvalid = False
    if "tRNA" in s:
        isvalid = False
    # if "ribosomal protein" in s:
    #     isvalid = False
    # if "rRNA" in s:
    #     isvalid = False
    # if "Ribosomal RNA" in s:
    #     isvalid = False
    return(isvalid)

def main():
    args = sys.argv
    if len(args) != 3:
        sys.exit("Insufficient arguments. Please specify path to genome")
    cdsfile = args[1]
    outdir = args[2]
    collect = []
    primer_counter = 0
    for record in SeqIO.parse(cdsfile,"fasta"):
        if is_valid(record.description):
            print(record.description)
            isvalid = False
            seq= get_primers(record)
            collect.extend(seq)
            primer_counter +=1
        if primer_counter > NUM_PRIMERS_TO_GENERATE:
            break
    df = pd.DataFrame(collect)
    genome = cdsfile.split("/")[-1].split(".")[0]
    df = df.assign(genome = genome)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print(df.shape)
    df = df.sort_values(by="amplicon_size")
    df.to_csv(f"{outdir}/{genome}.csv",index=False)
    with open(f"{outdir}/{genome}_R1.fasta","w") as outfile:
        for i, row in df.iterrows():
            outfile.write(f">{row.genome}:{row.gene}:{row.primer_id}:f\n{row.forward_primer}\n")
    with open(f"{outdir}/{genome}_R2.fasta","w") as outfile:
        for i, row in df.iterrows():
            outfile.write(f">{row.genome}:{row.gene}:{row.primer_id}:r\n{row.reverse_primer}\n")

if __name__ == '__main__':
    main()
