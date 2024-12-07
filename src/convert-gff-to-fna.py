import sys
from Bio import SeqIO, Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
    
def convert_to_fasta(genome, gff, output):
    Genome = next(SeqIO.parse(genome, "fasta"))
    collect = []
    GFF = pd.read_csv(gff, comment="#",
                      sep="\t",
                      names=["seqid","source","type","start","end", "score","strand","misc","annotation"])
    for i, feature in GFF.iterrows():
        if feature.strand == "+":
            strand = 1
        else:
            strand = -1
        print(feature.start, feature.end)
        f = SeqFeature(FeatureLocation(feature.start, feature.end, strand=strand),
                       type=feature["type"])
        seqid = ""
        if "=" in feature.annotation:
            seqid = feature.annotation.split(";")[0].split("=")[1]
        else:
            seqid = feature.annotation.split(";")[0].split(" ")[1]
        sr = SeqRecord(seq=f.extract(Genome).seq, features=[f], 
                       id=seqid, 
                       description=Genome.description)
        collect.append(sr)
    SeqIO.write(collect, output, "fasta")

def main():
    argv = sys.argv

    if len(argv) == 1:
        print("Missing arguments: Genome name, GFF file, output file name")
        sys.exit()

    genome = argv[1]
    gff = argv[2]
    output = argv[3]

    print(f"Genome: {genome}")
    convert_to_fasta(genome, gff, output)          


if __name__ == '__main__':
    main()
