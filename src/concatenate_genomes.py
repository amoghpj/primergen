import sys
from Bio import SeqIO, Seq
import os
def clean_fname(f):
    if "/" in f:
        f = f.split("/")[-1]
    if "." in f:
        f = f.split(".")[0]
    return(f)
    
def concatenate_genomes(flist):
    collect = []
    for f in flist:
        for record in SeqIO.parse(f,"fasta"):
            record.id = clean_fname(f) + "-" + record.id.replace(" ","_")
            collect.append(record)
    return(collect)

def main():
    argv = sys.argv
    if len(argv) == 1:
        print("Missing arguments: Genome name, and rest of genomes")
    genome = argv[1]
    output = argv[2]
    rest = argv[3:]
    genome = clean_fname(genome)
    print(f"Genome: {genome}")
    print(f"Rest: {rest}")
    concat = concatenate_genomes(rest)
    if not os.path.exists("concatenated/"):
        os.makedirs("concatenated/")
    
    SeqIO.write(concat, output, "fasta")

if __name__ == '__main__':
    main()
