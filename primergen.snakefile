import os
import glob
from pathlib import Path
import pandas as pd

run_path = Path(config["run_path"])

#### inputs
genome_dir = Path("genome/")

### outputs
annotation_dir = Path("cds/")
concat_dir = Path("concatenated/")
proposed_dir = Path("proposed/")

self_aligned_dir = Path("aligned/")
nonself_aligned_dir = Path("nonselfaligned/")
validate_nonself_primers_dir = Path("validate_nonself/")
self_index_dir = Path("genomeidx/")
nonself_index_dir = Path("concatenatedidx/")
target_dir = Path("target/")
output_dir = Path("design")

TARGETFILES = [f.split(".")[-2].split("/")[-1] for f in glob.glob(str(run_path / genome_dir) + "/*.fasta")]
EXCLUDEDGENOME = [ run_path / concat_dir / Path(f"{f}_excluded.fasta") for f in TARGETFILES]
ALLGENOMES = list(glob.glob(str(genome_dir) + "/*.fasta"))
ALLGENOMES.extend(list(glob.glob(str(concat_dir) + "/*.fasta")))

# print(TARGETFILES)

rule all:
    input:
        run_path / "top5_per_genome.csv"

rule aggregate:
    input:
        expand(run_path / output_dir / ("{genome}_primers.csv"), genome=TARGETFILES)
    output:
        run_path / "top5_per_genome.csv"
    run:
        dflist = []
        for f in input:
            dflist.append(pd.read_csv(f).head(5))
        df = pd.concat(dflist).reset_index(drop=True)
        df.to_csv(output, index=False)

rule find_cds:
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
        cds= run_path / annotation_dir / ("{genome}.ffn")
    conda:
        "/n/groups/springer/amogh/conda/bakta/"
    params:
        dbname="./db",
        outdir = run_path / annotation_dir
    shell:
        """mkdir -p {params.outdir}
bakta --db {params.dbname} {input.genome} -o {params.outdir} --skip-crispr --skip-plot --skip-pseudo --skip-trna --skip-ncrna --skip-rrna"""

rule propose_candidates:
    """
    Use primer3 to propose candidates for each genome
    """
    resources:
        runtime="10m",
        mem_mb="4000",
        partition="short",
    input:
        target= run_path / annotation_dir / ("{genome}.ffn")
    conda:
        "/n/groups/springer/amogh/conda/primer3/"
    output:
        proposed=run_path / proposed_dir /"{genome}.csv",
        R1 = run_path / proposed_dir /"{genome}_R1.fasta",
        R2 = run_path / proposed_dir /"{genome}_R2.fasta"
    params:
        outdir = os.path.join(run_path , proposed_dir)
    log:
        os.path.join(run_path, "logs", "{genome}_proposed.log")
        #lambda w : os.path.join(run_path ,"logs", w.genome + "_proposed.log")
    shell:
       """mkdir -p {params}/
python src/propose_primers.py {input.target} {params.outdir} >> {log}
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
bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2} -S {output} &> {log}
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
       run_path / "logs/{genome}_nonselfalign.log"
    shell:
       """mkdir -p {params.outputdir}/
bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2} -S {output} &> {log}
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
