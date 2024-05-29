import sys
from Bio import SeqIO, Seq
import os
import primer3 as p3
import pandas as pd

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
    # print(primers.keys())
    # sys.exit()
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
                               "amplicon_size":abs(startr + endr - (startf)),
                               "product_tm":primers[f"PRIMER_PAIR_{i}_PRODUCT_TM"],
                               "forward_primer":primers[f"PRIMER_LEFT_{i}_SEQUENCE"],
                               "reverse_primer":primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]})
        except:
            ok = False
        i += 1

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
    if "prfa: peptide chain release factor 1" in s:
        print(f"removing {s}")
        isvalid = False
    if "prfb: peptide chain release factor 2" in s:
        print(f"removing {s}")
        isvalid = False
    if "s20: ribosomal protein bS20" in s:
        print(f"removing {s}")
        isvalid = False
    if "tigr00043: rRNA maturation RNase YbeY" in s:
        print(f"removing {s}")
        isvalid = False
    if "tigr00054: RIP metalloprotease RseP" in s:
        print(f"removing {s}")
        isvalid = False
    if "l17: ribosomal protein bL17" in s:
        print(f"removing {s}")
        isvalid = False
    if "l21: ribosomal protein bL21" in s:
        print(f"removing {s}")
        isvalid = False
    if "ftsy: signal recognition particle-docking protein FtsY" in s:
        print(f"removing {s}")
        isvalid = False
    if "ftsz: cell division protein FtsZ" in s:
        print(f"removing {s}")
        isvalid = False
    if "rbfa: ribosome-binding factor A" in s:
        print(f"removing {s}")
        isvalid = False
    if "ribf: riboflavin biosynthesis protein RibF" in s:
        print(f"removing {s}")
        isvalid = False
    if "ruva: Holliday junction DNA helicase RuvA" in s:
        print(f"removing {s}")
        isvalid = False
    if "smpb: SsrA-binding protein" in s:
        print(f"removing {s}")
        isvalid = False
        "tif: tRNA (guanine(37)-N(1))-methyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "rsfs_iojap_ybeB: ribosome silencing factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "tigr00092: GTP-binding protein YchF" in s:
        print(f"removing {s}")
        isvalid = False
    if "tigr00095: 16S rRNA (guanine(966)-N(2))-methyltransferase RsmD" in s:
        print(f"removing {s}")
        isvalid = False
    if "tig: trigger factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "tsf: translation elongation factor Ts" in s:
        print(f"removing {s}")
        isvalid = False
    if "rsmg_gidB: 16S rRNA (guanine(527)-N(7))-methyltransferase RsmG" in s:
        print(f"removing {s}")
        isvalid = False
    if "l9: ribosomal protein bL9" in s:
        print(f"removing {s}")
        isvalid = False
    if "s6: ribosomal protein bS6" in s:
        print(f"removing {s}")
        isvalid = False
    if "infc: translation initiation factor IF-3" in s:
        print(f"removing {s}")
        isvalid = False
    if "rrna_methyl_3: RNA methyltransferase, TrmH family, group 3" in s:
        print(f"removing {s}")
        isvalid = False
    if "uvrc: excinuclease ABC subunit C" in s:
        print(f"removing {s}")
        isvalid = False
        "rife_H_YqgF: putative transcription antitermination factor YqgF" in s:
        print(f"removing {s}")
        isvalid = False
    if "pyrg: CTP synthase" in s:
        print(f"removing {s}")
        isvalid = False
    if "alas: alanine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "dnaa: chromosomal replication initiator protein DnaA" in s:
        print(f"removing {s}")
        isvalid = False
    if "clpx: ATP-dependent Clp protease, ATP-binding subunit ClpX" in s:
        print(f"removing {s}")
        isvalid = False
    if "iles: isoleucine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "leus_bact: leucine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "metg: methionine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "sers: serine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "sms: DNA repair protein RadA" in s:
        print(f"removing {s}")
        isvalid = False
    if "trmu: tRNA (5-methylaminomethyl-2-thiouridylate)-methyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "trub: tRNA pseudouridine(55) synthase" in s:
        print(f"removing {s}")
        isvalid = False
        "cif: cysteine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
        "era: GTifinding protein Era" in s:
        print(f"removing {s}")
        isvalid = False
        "hiss: histidiif-tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "mray: phospho-N-acetylmuramoyl-pentapeptide-transferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "args: arginine--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "asps_bact: aspartate--tRNA ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "fmt: methionyl-tRNA formyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "phes: phenylalanine--tRNA ligase, alpha subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "phet_bact: phenylalanine--tRNA ligase, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "if-2: translation initiation factor IF-2" in s:
        print(f"removing {s}")
        isvalid = False
    if "frr: ribosome recycling factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "hemn_rel: putative oxygen-independent coproporphyrinogen III oxidase" in s:
        print(f"removing {s}")
        isvalid = False
    if "mfd: transcription-repair coupling factor" in s:
        print(f"removing {s}")
        isvalid = False
    if "pola: DNA polymerase I" in s:
        print(f"removing {s}")
        isvalid = False
        "rif: recombination protein RecR" in s:
        print(f"removing {s}")
        isvalid = False
        "uvrb: eifnuclease ABC subunit B" in s:
        print(f"removing {s}")
        isvalid = False
        "recn: DNA repif protein RecN" in s:
        print(f"removing {s}")
        isvalid = False
    if "ruvb: Holliday junction DNA helicase RuvB" in s:
        print(f"removing {s}")
        isvalid = False
    if "recg: ATP-dependent DNA helicase RecG" in s:
        print(f"removing {s}")
        isvalid = False
    if "dnan: DNA polymerase III, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpsa: ribosomal protein bS1" in s:
        print(f"removing {s}")
        isvalid = False
    if "ksga: ribosomal RNA small subunit methyltransferase A" in s:
        print(f"removing {s}")
        isvalid = False
    if "secg: preprotein translocase, SecG subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "nusg: transcription termination/antitermination factor NusG" in s:
        print(f"removing {s}")
        isvalid = False
    if "purb: adenylosuccinate lyase" in s:
        print(f"removing {s}")
        isvalid = False
    if "ffh: signal recognition particle protein" in s:
        print(f"removing {s}")
        isvalid = False
    if "seca: preprotein translocase, SecA subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "sece_bact: preprotein translocase, SecE subunit" in s:
        print(f"removing {s}")
        isvalid = False
        "3if01s007: preprotein translocase, SecY subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpsc_bact: ribosomal protein uS3" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpsb_bact: ribosomal protein uS2" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpsd_bact: ribosomal protein uS4" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpse_bact: ribosomal protein uS5" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpsg_bact: ribosomal protein uS7" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplt_bact: ribosomal protein bL20" in s:
        print(f"removing {s}")
        isvalid = False
    if "atpd: ATP synthase F1, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplv_bact: ribosomal protein uL22" in s:
        print(f"removing {s}")
        isvalid = False
    if "gyrb: DNA gyrase, B subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "gyra: DNA gyrase, A subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplm_bact: ribosomal protein uL13" in s:
        print(f"removing {s}")
        isvalid = False
        "rif_bact: ribosomal protein uL15" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplx_bact: ribosomal protein uL24" in s:
        print(f"removing {s}")
        isvalid = False
    if "murc: UDP-N-acetylmuramate--L-alanine ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "murd: UDP-N-acetylmuramoylalanine--D-glutamate ligase" in s:
        print(f"removing {s}")
        isvalid = False
    if "hola: DNA polymerase III, delta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "atpsyn_F1gamma: ATP synthase F1, gamma subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplp_bact: ribosomal protein uL16" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpla_bact: ribosomal protein uL1" in s:
        print(f"removing {s}")
        isvalid = False
    if "rplb_bact: ribosomal protein uL2" in s:
        print(f"removing {s}")
        isvalid = False
    if "imp_dehydrog: inosine-5'-monophosphate dehydrogenase" in s:
        print(f"removing {s}")
        isvalid = False
    if "dnag: DNA primase" in s:
        print(f"removing {s}")
        isvalid = False
    if "lepa: elongation factor 4" in s:
        print(f"removing {s}")
        isvalid = False
        "tif_BipA: GTP-binding protein TypA/BipA" in s:
        print(f"removing {s}")
        isvalid = False
    if "coad_prev_kdtB: pantetheine-phosphate adenylyltransferase" in s:
        print(f"removing {s}")
        isvalid = False
    if "l11_bact: ribosomal protein uL11" in s:
        print(f"removing {s}")
        isvalid = False
    if "nusb: transcription antitermination factor NusB" in s:
        print(f"removing {s}")
        isvalid = False
    if "nusa: transcription termination factor NusA" in s:
        print(f"removing {s}")
        isvalid = False
    if "tigrfam_recA: protein RecA" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpob: DNA-directed RNA polymerase, beta subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpoa: DNA-directed RNA polymerase, alpha subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "pyrh_bact: UMP kinase" in s:
        print(f"removing {s}")
        isvalid = False
    if "rnaseiii: ribonuclease III" in s:
        print(f"removing {s}")
        isvalid = False
    if "16s_RimM: 16S rRNA processing protein RimM" in s:
        print(f"removing {s}")
        isvalid = False
    if "prok_dnaK: chaperone protein DnaK" in s:
        print(f"removing {s}")
        isvalid = False
    if "rif_TIGR: DNA-directed RNA polymerase, beta' subunit" in s:
        print(f"removing {s}")
        isvalid = False
    if "dnax_ntif: DNA polymerase III, subunit gamma and tau" in s:
        print(f"removing {s}")
        isvalid = False
    if "lysidine_TilS_N: tRNA(Ile)-lysidine synthetase" in s:
        print(f"removing {s}")
        isvalid = False
    if "obg_CgtA: Obg family GTPase CgtA" in s:
        print(f"removing {s}")
        isvalid = False
    if "guanyl_kin: guanylate kinase" in s:
        print(f"removing {s}")
        isvalid = False
    if "gtpase_EngA: ribosome-associated GTPase EngA" in s:
        print(f"removing {s}")
        isvalid = False
    if "l3_bact: 50S ribosomal protein uL3" in s:
        print(f"removing {s}")
        isvalid = False
    if "us11_bact: ribosomal protein uS11" in s:
        print(f"removing {s}")
        isvalid = False
    if "l6_bact: ribosomal protein uL6" in s:
        print(f"removing {s}")
        isvalid = False
    if "t6a_TsaD_YgjD: tRNA threonylcarbamoyl adenosine modification protein TsaD" in s:
        print(f"removing {s}")
        isvalid = False
    if "t6a_YeaZ: tRNA threonylcarbamoyl adenosine modification protein YeaZ" in s:
        print(f"removing {s}")
        isvalid = False
    if "rpld_bact: 50S ribosomal protein uL4" in s:
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
    if len(args) != 2:
        sys.exit("Insufficient arguments. Please specify path to genome")
    cdsfile = args[1]
    collect = []
    counter = 0
    for record in SeqIO.parse(cdsfile,"fasta"):
        if is_valid(record.description):
            print(record.description)
            isvalid = False
            seq= get_primers(record)
            collect.extend(seq)
            counter +=1
        if counter > 100:
            break
    df = pd.DataFrame(collect)
    genome = cdsfile.split("/")[-1].split(".")[0]
    df = df.assign(genome = genome)
    if not os.path.exists("proposed/"):
        os.makedirs("proposed/")
    print(df.shape)
    df = df.sort_values(by="amplicon_size")
    df.to_csv(f"proposed/{genome}.csv",index=False)
    with open(f"proposed/{genome}_R1.fasta","w") as outfile:
        for i, row in df.iterrows():
            outfile.write(f">{row.genome}:{row.gene}:{row.primer_id}:f\n{row.forward_primer}\n")
    with open(f"proposed/{genome}_R2.fasta","w") as outfile:
        for i, row in df.iterrows():
            outfile.write(f">{row.genome}:{row.gene}:{row.primer_id}:r\n{row.reverse_primer}\n")

if __name__ == '__main__':
    main()
