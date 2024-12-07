import os
import glob
from pathlib import Path
import pandas as pd

run_path = Path(config["run_path"])

#### inputs
genome_dir = Path("genome/")
### OPTIONAL
if config.get("strain_type", None) is not None:
    straintypedf = pd.read_csv(Path(config["strain_type"]))
    STRAINDICT = {row.strain: row["type"] for i,row in straintypedf.iterrows()}
    print(STRAINDICT)

### outputs
annotation_dir = Path("cds/")
concat_dir = Path("concatenated/")
proposed_dir = Path("proposed/")

self_aligned_dir = Path("aligned_long/")
nonself_aligned_dir = Path("nonselfaligned_long/")
validate_nonself_primers_dir = Path("validate_nonself_long/")
self_index_dir = Path("genomeidx/")
nonself_index_dir = Path("concatenatedidx/")
target_dir = Path("target/")
output_dir = Path("design_long")

TARGETFILES = [f.split(".")[-2].split("/")[-1] for f in glob.glob(str(run_path / genome_dir) + "/*.fasta")]
EXCLUDEDGENOME = [ run_path / concat_dir / Path(f"{f}_excluded.fasta") for f in TARGETFILES]
ALLGENOMES = list(glob.glob(str(genome_dir) + "/*.fasta"))
ALLGENOMES.extend(list(glob.glob(str(concat_dir) + "/*.fasta")))

# print(TARGETFILES)

rule all:
    input:
        run_path / "top10_per_genome_long.csv"

rule aggregate:
    resources:
        runtime="10m",
        mem_mb="500",
        partition="short",
    input:
        expand(run_path / output_dir / ("{genome}_primers.csv"), genome=TARGETFILES)
    output:
        run_path / "top10_per_genome_long.csv"
    run:
        dflist = []
        noprimers = []
        for f in input:
            try: 
                df = pd.read_csv(f)
                dflist.append(df.groupby("gene").first().reset_index().head(10))
            except:
                noprimers.append(f)
        df = pd.concat(dflist).reset_index(drop=True)
        df.to_csv(run_path / "top10_per_genome_long.csv", index=False)
        with open("no_primers_found","w") as outfile:
            outfile.write("\n".join(noprimers))

rule find_cds_prokaryote:
    """
    If CDS file for a genome is absent, use $annotation_tool to
    predict coding sequences and annotate them.
    """
    resources:
        runtime="3h",
        mem_mb="16000",
        cpus_per_task=4,
        partition="short",
    input:
        genome= run_path / genome_dir/ ("{genome}.fasta")
    output:
        cds= run_path / annotation_dir / ("{genome}_prokaryote.ffn")
    conda:
        "/n/groups/springer/amogh/conda/bakta/"
    params:
        dbname="./db",
        outdir = run_path / annotation_dir,
        rundir = lambda wc : run_path / Path("prokaryotic_genomes") / wc.genome ,
        outffn = "{genome}.ffn",
    shell:
        """mkdir -p {params.outdir}
mkdir -p {params.rundir}
bakta --db {params.dbname} {input.genome} -o {params.rundir} --skip-crispr --skip-plot --skip-pseudo --skip-trna --skip-ncrna --skip-rrna
cp {params.rundir}/{params.outffn} {output.cds}
"""


rule find_cds_eukaryote:
    """
    If CDS file for a genome is absent, use $annotation_tool to
    predict coding sequences and annotate them.
    """
    resources:
        runtime="4h",
        mem_mb="16000",
        cpus_per_task=4,
        partition="short",
    input:
        genome= run_path / genome_dir/ ("{genome}.fasta")
    output:
        cds= run_path / annotation_dir / ("{genome}_eukaryote.ffn")
    params:
        outdir = run_path / annotation_dir,
        rundir = lambda wc : run_path / Path("eukaryotic_genomes") / Path(wc.genome),
        inputgenome = lambda wc : (wc.genome + ".fasta")
    shell:
        """mkdir -p {params.outdir}/
mkdir -p {params.rundir}/
cp {input.genome} {params.rundir}
cd {params.rundir}
perl /n/groups/springer/amogh/src/genemarks/gmes_linux_64/gmes_petap.pl --sequence {params.inputgenome} --ES
cd ../../../
python src/convert-gff-to-fna.py '{input.genome}' '{params.rundir}/genemark.gtf' '{output.cds}'
"""

rule find_cds_virus:
    """
    If CDS file for a genome is absent, use $annotation_tool to
    predict coding sequences and annotate them.
    """
    resources:
        runtime="15m",
        mem_mb="1G",
        cpus_per_task=1,
        partition="short",
    conda:
        "/n/groups/springer/amogh/conda/seq/"
    input:
        genome= run_path / genome_dir/ ("{genome}.fasta")
    output:
        cds= run_path / annotation_dir / ("{genome}_virus.ffn")
    params:
        outdir = run_path / annotation_dir,
        rundir = lambda wc : run_path / Path("virus_genomes") / Path(wc.genome),
        outgff = "{genome}.gff"
    shell:
        """mkdir -p {params.outdir}
mkdir -p {params.rundir}
prodigal -i {input.genome} -f gff -o {params.rundir}/{params.outgff} -p meta
python src/convert-gff-to-fna.py {input.genome} {params.rundir}/{params.outgff} {output.cds}
"""

rule propose_candidates:
    """
    Use primer3 to propose candidates for each genome
    """
    resources:
        runtime="20m",
        mem_mb="4000",
        partition="short",
    input:
        target = lambda wc : run_path / annotation_dir / (wc.genome + "_" + STRAINDICT[wc.genome] + ".ffn")
    conda:
        "/n/groups/springer/amogh/conda/primer3/"
    output:
        proposed=run_path / proposed_dir /"{genome}.csv",
        R1 = run_path / proposed_dir /"{genome}_R1.fasta",
        R2 = run_path / proposed_dir /"{genome}_R2.fasta"
    params:
        outdir = os.path.join(run_path , proposed_dir),
        prefix = "{genome}"
    log:
        os.path.join(run_path, "logs", "{genome}_proposed.log")
        #lambda w : os.path.join(run_path ,"logs", w.genome + "_proposed.log")
    shell:
       """mkdir -p {params}/
python src/propose_primers.py {input.target} {params.outdir} {params.prefix} >> {log}
"""

rule concatenate_genomes:
    """For every genome in GENOMELIST, exclude it, and make a
    concatenated genome FASTA file for every other genome. 
    """
    resources:
        runtime="5m",
        mem_mb="6000",
        partition="short",
    input:
        genome= run_path / genome_dir/ "{genome}.fasta"
    output:
        excluded= run_path / concat_dir / "{genome}_excluded.fasta"
    params:
        genome= lambda w : w.genome
    log:
        os.path.join(run_path, "logs", "{genome}_concatenate.log")
    run:
        excluded = []
        for f in list(TARGETFILES):
            if f != params.genome:
                excluded.append(f"{str(run_path)}/genome/{f}.fasta")
        shell(f"python src/concatenate_genomes.py {input.genome} {output.excluded} {' '.join(excluded)} > {log}")
    
rule index_genomes:
    resources:
        runtime="30m",
        mem_mb="10000",
        partition="short",
    conda:
        "/n/groups/springer/amogh/conda/bowtie2/"
    input:
        genome= run_path / genome_dir / "{genome}.fasta",
        excludedgenome= run_path / concat_dir / "{genome}_excluded.fasta"
    output:
        genomeidx= run_path / self_index_dir / "{genome}.1.bt2",
        exclgenomeidx= run_path / nonself_index_dir / "{genome}_excluded.1.bt2",
    params:
        genome= lambda w: w.genome,
        genomeidx= lambda w : str(run_path / self_index_dir) + "/" + w.genome,
        exclgenomeidx= lambda w : str(run_path / nonself_index_dir) + "/" + w.genome + "_excluded",
        selfidxdir = run_path / self_index_dir,
        nonselfidxdir = run_path / nonself_index_dir,
    log:
        os.path.join(run_path, "logs", "{genome}_index.log")
    shell:
        """mkdir -p {params.selfidxdir}/
mkdir -p {params.nonselfidxdir}/
bowtie2-build {input.genome} {params.genomeidx} 1> {log}
bowtie2-build {input.excludedgenome} {params.exclgenomeidx} 1>> {log}
"""

rule align_self:
    resources:
        runtime="10m",
        mem_mb="2000",
        partition="short",
    input:
        R1= run_path / proposed_dir / "{genome}_R1.fasta",
        R2=run_path / proposed_dir / "{genome}_R2.fasta",
        index= run_path / self_index_dir / "{genome}.1.bt2"
    output:
        samfile= run_path / self_aligned_dir/ "{genome}.sam",
    params:
        indexref= lambda w : str(run_path / self_index_dir) + "/" +  w.genome,
        outputdir = run_path / self_aligned_dir
    conda:
        "/n/groups/springer/amogh/conda/bowtie2/"
    log:
       run_path / "logs/{genome}_align.log"
    shell:
       """mkdir -p {params.outputdir}/
bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2} -I 90 -X 2000 -S {output} &> {log}
"""

rule align_nonself:
    resources:
        runtime="10m",
        mem_mb="2000",
        partition="short",
    input:
        R1= run_path / proposed_dir / "{genome}_R1.fasta",
        R2=run_path / proposed_dir / "{genome}_R2.fasta",
        index= run_path / nonself_index_dir / "{genome}_excluded.1.bt2"
    output:
        samfile= run_path / nonself_aligned_dir/ "{genome}.sam",
    conda:
        "/n/groups/springer/amogh/conda/bowtie2/"
    params:
        indexref= lambda w : str(run_path / nonself_index_dir) + "/" +  w.genome + "_excluded",
        outputdir = run_path / nonself_aligned_dir
    log:
       run_path / "logs/{genome}_nonselfalign_long.log"
    shell:
       """mkdir -p {params.outputdir}/
bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2}  -I 0 -X 5000 -S {output} &> {log}
"""

rule filter_candidates:
    resources:
        runtime="5m",
        mem_mb="500",
        partition="short",
    input:
        proposed=run_path / proposed_dir / "{genome}.csv",
        aligntocheck= run_path / nonself_aligned_dir / "{genome}.sam",
        sanitycheck= run_path / self_aligned_dir /"{genome}.sam"
    output:
        run_path / output_dir/ "{genome}_primers.csv"
    params:
        output = run_path / output_dir
    log:
       run_path / "logs/{genome}_filter.log"
    shell:
        """mkdir -p {params.output}/
python src/filter_primers.py {input.proposed} {input.aligntocheck} {output} 1> {log}"""
