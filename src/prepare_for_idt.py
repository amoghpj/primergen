import os
import sys
import time
import datetime
import pandas as pd
import primer3 as p3
from Bio import SeqIO, Seq

fluorquench=[["FAM", "5' 6-FAM", "/56-FAM/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["TET", "5' TET", "/5TET/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["Yakima", "5' Yakima Yellow®", "/5YakYel/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["HEX", "5' HEX", "/5HEX/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["JOE", "5' JOE™ (NHS Ester)", "/56-JOEN/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["Cy3", "5' Cy®3", "/5Cy3/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["TexasRedX", "5' Texas Red®-X (NHS Ester)", "/5TexRd-XN/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["Cy5", "5' Cy®5", "/5Cy5/", "TAO - 3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["MAX", "5' MAX™ (NHS Ester)", "/5MAXN/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["Tye563", "5' TYE™ 563", "/5TYE563/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["TAMRA", "5' TAMRA (NHS Ester)", "/56-TAMN/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["ROX", "5' ROX (NHS Ester)", "/56-ROXN/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["TEX615", "5' TEX 615™", "/5TEX615/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["TYE625", "5' TYE™ 665", "/5TYE665/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["ATTO425", "5' ATTO™ 425 (NHS Ester)", "/5ATTO425N/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"],
             ["ATTO550", "5' ATTO™ 550 (NHS Ester)", "/5ATTO550N/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["ATTO590", "5' ATTO™ 590 (NHS Ester)", "/5ATTO590N/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["ATTO647N", "5' ATTO™ 647N (NHS Ester)", "/5ATTO647NN/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["ATTO700", "5' ATTO™ 700 (NHS Ester)", "/5ATTO700N/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["Cy5.5", "5' Cy®5.5", "/5Cy55/", "3' Iowa Black™ RQ-Sp", "/3IAbRQSp/"],
             ["SUN", "5' SUN™", "/5SUN/", "ZEN - 3' Iowa Black™ FQ", "/3IABkFQ/"]]

fluordf = pd.DataFrame(fluorquench,
                       columns=["5' Dye","5' Dye Name",
                                "5' Code", "3' Quencher",
                                "3' Code"])

def transform_probe_sequence(seq, fluor):
    """
    Add 5' dye, 3' quencher, and optionally a second quencher (Zen or Tao)
    """
    fluordetails = fluordf[fluordf["5' Dye"] == fluor]
    code_5,code_3 = fluordetails["5' Code"].values[0],fluordetails["3' Code"].values[0]
    secondquencher = ""
    quencher = fluordetails["3' Quencher"].values[0].split("-")[0].strip()
    if quencher in ["ZEN", "TAO"]:
        secondquencher = quencher
        seq = seq[:9] + "/" + secondquencher + "/" + seq[9:]
        modifiedseq = code_5 + seq + code_3
    return(modifiedseq)

def to_excel(df, outpath):
    writer = pd.ExcelWriter(outpath, engine='xlsxwriter')
    df.to_excel(writer, index=False, sheet_name='Sheet1')
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    writer.close()

def main():
    args = sys.argv
    if len(args) != 4:
        sys.exit("Insufficient arguments. Please specify input path, output path, and optionally the straintype file.")

    inpath = args[1]
    runname = args[2]
    prefixpath = args[3]
    exporttemplate = pd.read_excel("./data/templates/IDT_primer_template.xlsx")
    inputdf = pd.read_csv(inpath)\
             .groupby("genome").head(1)\
                               .reset_index()\
                               [["genome","gene","primer_id",
                                 "forward_tm","reverse_tm","product_tm",
                                 "internal_sequence",
                                 "forward_primer","reverse_primer"]]
    straindf = pd.read_csv(prefixpath)
    if "prefix" in straindf.columns:
        inputdf = inputdf.merge(straindf, left_on="genome",right_on="strain")
    else:
        inputdf = inputdf.assign(prefix=[v.split("_")[0][:5] + "_" + v.split("_")[1][:5]\
                                 for v in inputdf.genome.values])

    #inputdf = inputdf.assign(prefix_abbrev = inputdf.prefix.str.cat(inputdf.abbrev, sep=":"))
    # print(inputdf.columns)
    # print(inputdf.prefix)
    for i, row in inputdf.iterrows():
        for oligo, shorthand in zip(["forward_primer",
                                     "reverse_primer",
                                     "internal_sequence"],
                                    ["F","R","P"]):
            if "primer" in oligo:
                exporttemplate = pd.concat([exporttemplate,
                                            pd.DataFrame([{"Name":f"{row.prefix}_{shorthand}",
                                                           "Sequence":row[oligo],
                                                           "Scale":"250nm",
                                                           "Purification":"STD"}])])
            else:
                """
                Assume the default fluor is FAM. 
                TODO allow user to specify fluor in config file.
                """
                fluor = "FAM"
                exporttemplate = pd.concat([exporttemplate,
                                            pd.DataFrame([{"Name":f"{row.prefix}_{shorthand}",
                                                           "Sequence": transform_probe_sequence(row[oligo], fluor),
                                                           "Scale":"250nm",
                                                           "Purification":"HPLC"}])])
    export_path = f"{runname}/{runname}_idt.xlsx"
    to_excel(exporttemplate, export_path)

if __name__ == '__main__':
    main()
