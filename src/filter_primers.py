import sys
import os 
import pandas as pd

def main():
    argv = sys.argv
    if len(argv) == 1:
        print("Missing primer list")
        sys.exit()
    sam = argv[2]
    proposed = argv[1]
    output = argv[3]
    samdf = pd.read_csv(sam,
                        sep="\t",
                        comment="@",
                        names=["QNAME","FLAG","RNAME","POS","MAPQ","CIGAR",
                               "MRNM","MPOS","ISIZE","SEQ","QUAL1","QUAL2","QUAL3","QUAL4","QUAL5","QUAL6","QUAL7","QUAL8","QUAL9","QUAL10","QUAL11","QUAL12",],
                        engine="python")
    samdf = samdf[samdf.RNAME == "*"]
    primers = samdf.QNAME.str.split(":", expand=True).rename({0:"genome",1:"gene",2:"primer_id",3:"orientation"}, axis="columns")
    if primers.shape[0] == 0:
        pd.DataFrame({}).to_csv(output, index=False)
    else:
        primers["primer_id"] = primers.primer_id.astype(int)
        primers = primers[["genome","gene","primer_id"]]
        primers = primers.drop_duplicates().reset_index(drop=True)
        print(primers)
        proposeddf = pd.read_csv(proposed)
        outdf = primers.merge(proposeddf, on=["genome","gene","primer_id"])
        print("Input count:", proposeddf.shape[0])
        print("Filtered count:", outdf.shape[0])
        outdf.to_csv(output, index=False)


if __name__ == "__main__":
    main()
