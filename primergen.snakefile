import os
import glob
from pathlib import Path

run_path = Path(config["run_path"])
#### inputs

genome_dir = Path("genome/")
annotation_dir = Path("cds/")

### outputs
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

print(TARGETFILES)

rule all:
    input:
        expand(run_path / output_dir / ("{genome}_primers.csv"), genome=TARGETFILES)

rule find_cds:
    """
    If CDS file for a genome is absent, use $annotation_tool to
    predict coding sequences and annotate them.
    """
    input:
        genome= run_path / genome_dir/ ("{genome}.fasta")
    output:
        cds= run_path / annotation_dir / ("{genome}.ffn")
    conda:
        "/n/groups/springer/amogh/conda/bakta/"
    params:
        dbname="./BAKTA_DB",
        outdir = run_path / annotation_dir
    shell:
        """mkdir -p {params.outdir}
bakta --db {params.dbname} {input.genome} {output.cds}"""

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
    output:
        proposed=run_path / proposed_dir /"{genome}.csv",
        R1 = run_path / proposed_dir /"{genome}_R1.fasta",
        R2 = run_path / proposed_dir /"{genome}_R2.fasta"
    log:
        run_path / "logs/{genome}_proposed.log"
    run:
       shell("mkdir -p {{run_path}}/proposed/")
       shell(f"python src/propose_primers.py {input.target} > {log}")

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
        genome="{genome}"
    log:
       run_path / "logs/{genome}_concatenate.log"
    run:
        excluded = []
        for f in list(TARGETFILES):
            if f != params.genome:
                excluded.append(f"{str(run_path)}/genome/{f}.fasta")
        shell(f"python src/concatenate_genomes.py {input.genome} {' '.join(excluded)} > {log}")
    
rule index_genomes:
    resources:
        runtime="30m",
        mem_mb="10000",
        partition="short",
    input:
        genome= run_path / genome_dir / "{genome}.fasta",
        excludedgenome= run_path / concat_dir / "{genome}_excluded.fasta"
    output:
        genomeidx= run_path / self_index_dir / "{genome}.1.bt2",
        exclgenomeidx= run_path / nonself_index_dir / "{genome}_excluded.1.bt2",
    params:
        genome="{genome}",
        genomeidx=run_path / self_index_dir / "{genome}",
        exclgenomeidx=run_path / nonself_index_dir / "{genome}_excluded",
    log:
       run_path / ("logs/{genome}_index.log")
    run:
        shell("mkdir -p {params.genomeidx}/")
        shell("mkdir -p {params.exclgenomeidx}")
        shell(f"bowtie2-build {input.genome} {params.genomeidx} 1> {params.genome}.log")
        shell(f"bowtie2-build {input.excludedgenome} {params.exclgenomeidx} 1> {params.genome}.log")

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
        indexref= run_path / self_index_dir / "{genome}",
        outputdir = run_path / self_aligned_dir
    log:
       run_path / "logs/{genome}_align.log"
    run:
       shell("mkdir -p {params.outputdir}/")
       shell(f"bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2} -S {output} &> {log}")

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
    params:
        indexref= run_path / nonself_index_dir / "{genome}_excluded",
        outputdir = run_path / nonself_aligned_dir
    log:
       run_path / "logs/{genome}_nonselfalign.log"
    run:
       shell("mkdir -p {params.outputdir}/")
       shell(f"bowtie2 -x {params.indexref} -f -1 {input.R1} -2 {input.R2} -S {output} &> {log}")


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
    run:
        shell("mkdir -p {params.output}/")
        shell(f"python src/filter_primers.py {input.proposed} {input.aligntocheck} {output} 1> {log}")