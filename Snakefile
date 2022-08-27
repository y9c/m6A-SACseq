import sys
from collections import defaultdict

workdir: config["workdir"]

bed_number = config["bed_number"]
src_dir = config["srcdir"]

data_dir = os.path.dirname(workflow.overwrite_configfiles[-1])
ref_dir = data_dir

if "samples" not in config:
    sys.exit("`samples` is not defined in config file!")

if "references" not in config:
    sys.exit("`references` is not defined in config file!")

REF = config["references"]
# print(workflow.basedir)




group2sample = defaultdict(list)
run2file = {}
sample_ids = []
group_ids = []
run_ids = []
sample2run = defaultdict(list)
read_ids = set()
# group: eg, WT, KD
for g, v in config["samples"].items():
    group_ids.append(g)
    # lib: input, treated
    for l, v2 in v.items():
        for r, files in v2.items():
            s = "-".join([g, l, r])
            sample_ids.append(s)
            group2sample[g].append(s)
            for i, r in enumerate(files, 1):
                run_ids.append(s + f"-run{i}")
                sample2run[s].append(s + f"-run{i}")
                run2file[s + f"-run{i}"] = r
                read_ids |= set(r.keys())
# make sure R1 and R2 are in the correct order
read_ids = sorted(list(read_ids))

rule all:
    input:
        # before combine
        "quality_control/report_falco_before.html",
        "quality_control/report_cutadapt.html",
        "quality_control/report_falco_after.html",
        "quality_control/report_star.html",
        "quality_control/report_falco_unmap.html",
        # after combine
        "quality_control/report_rnaseqc.html",
        expand("stat_reads/{sample}.tsv", sample=sample_ids),
        expand("stat_mapping/{sample}.tsv", sample=sample_ids),
        expand("stat_dedup/{sample}.tsv", sample=sample_ids),
        expand("spike_aligned/{sample}.tsv.gz", sample=sample_ids),
        expand(
            "count_depth_by_sample/{sample}_{reftype}.tsv.gz",
            sample=sample_ids,
            reftype=["genome"],
        ),
        "count_reads/genome_single.count",
        expand(
            "pileup_filtered_by_group/{group}_{reftype}.tsv.gz",
            group=group_ids,
            reftype=["genome"],
        ),
        # reftype=["rRNA", "smallRNA", "genome"],


## CHANGE: do not combine runs before mapping, so we do not need to re-run the mapping when add new runs


# pre-QC


rule falco_before:
    input:
        lambda wildcards: os.path.join(data_dir, run2file[wildcards.rn][wildcards.rd]),
    output:
        html="quality_control/falco_before/{rn}_{rd}/fastqc_report.html",
        text="quality_control/falco_before/{rn}_{rd}/fastqc_data.txt",
        summary="quality_control/falco_before/{rn}_{rd}/summary.txt",
    params:
        "quality_control/falco_before/{rn}_{rd}",
    shell:
        "falco -o {params} {input}"


rule report_falco_before:
    input:
        reports=expand(
            "quality_control/falco_before/{rn}_{rd}/fastqc_data.txt",
            rn=run_ids,
            rd=read_ids,
        ),
    output:
        "quality_control/report_falco_before.html",
    resources:
        mem="20G",
    shell:
        "multiqc -f -m fastqc -n {output} {input.reports}"


# modify seq


rule cutadapt:
    input:
        lambda wildcards: [
            os.path.join(data_dir, run2file[wildcards.rn][rd]) for rd in read_ids
        ],
    output:
        inter_1=temp("cut_adapter/{rn}_cut_R1.fq"),
        inter_2=temp("cut_adapter/{rn}_cut_R2.fq"),
        trimmed_1="cut_adapter/{rn}_cut_R1.fq.gz",
        trimmed_2="cut_adapter/{rn}_cut_R2.fq.gz",
        # untrimmed_1="cut_adapter/{rn}_untrimmed_R1.fq.gz",
        # untrimmed_2="cut_adapter/{rn}_untrimmed_R2.fq.gz",
        short_1="cut_adapter/{rn}_short_R1.fq.gz",
        short_2="cut_adapter/{rn}_short_R2.fq.gz",
        report1="cut_adapter/{rn}_step1.report",
        report2="cut_adapter/{rn}_step2.report",
    params:
        adapter3_r1=config["barcode"]["barcode3_r1"],
        adapter3_r2=config["barcode"]["barcode3_r2"],
        primerF=config["barcode"]["RT_primer_F"],
        primerR=config["barcode"]["RT_primer_R"],
    threads: 20
    resources:
        mem="30G",
    shell:
        """
        cutadapt -j {threads} \
          -U 11 \
          --rename='{{id}}_{{r2.cut_prefix}} {{comment}}' \
          --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 \
          -O 6 \
          --pair-filter=both \
          -a {params.adapter3_r1} -A {params.adapter3_r2} \
          -o {output.inter_1} -p {output.inter_2} \
          {input} >{output.report1} 
        cutadapt -j {threads} \
          -m 15 \
          -u -11 \
          -n 5 \
          -O 12 \
          -g {params.primerF} -a {params.primerR} \
          -G {params.primerF} -A {params.primerR} \
          --too-short-output={output.short_1} --too-short-paired-output={output.short_2} \
          -o {output.trimmed_1} -p {output.trimmed_2} \
          {output.inter_1} {output.inter_2} >{output.report2} 
        """


rule report_cutadapt:
    input:
        reports=expand("cut_adapter/{rn}_step1.report", rn=run_ids),
    output:
        "quality_control/report_cutadapt.html",
    shell:
        "multiqc -f -m cutadapt -n {output} {input.reports}"


# post-QC


rule falco_after:
    input:
        "cut_adapter/{rn}_cut_{rd}.fq.gz",
    output:
        html="quality_control/falco_after/{rn}_{rd}/fastqc_report.html",
        text="quality_control/falco_after/{rn}_{rd}/fastqc_data.txt",
        summary="quality_control/falco_after/{rn}_{rd}/summary.txt",
    params:
        "quality_control/falco_after/{rn}_{rd}",
    shell:
        "falco -o {params} {input}"


rule report_falco_after:
    input:
        reports=expand(
            "quality_control/falco_after/{rn}_{rd}/fastqc_data.txt",
            rn=run_ids,
            rd=read_ids,
        ),
    output:
        "quality_control/report_falco_after.html",
    resources:
        mem="20G",
    shell:
        "multiqc -f -m fastqc -n {output} {input.reports}"


## Mapping: 1, contamination


rule map_to_contamination_by_bowtie2:
    input:
        r1="cut_adapter/{rn}_cut_R1.fq.gz",
        r2="cut_adapter/{rn}_cut_R2.fq.gz",
    output:
        sam=temp("bowtie2_mapping/{rn}_contamination.sam"),
        un1=temp("bowtie2_mapping/{rn}_contamination.1.fq"),
        un2=temp("bowtie2_mapping/{rn}_contamination.2.fq"),
        report="bowtie2_mapping/{rn}_contamination.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: os.path.join(ref_dir, REF["contamination"]["bt2"]),
        un="bowtie2_mapping/{rn}_contamination.fq",
    threads: 24
    resources:
        mem="96G",
    shell:
        # The condition is not that stringent
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --no-unal --end-to-end --fast --un-conc {params.un} -x {params.ref_bowtie2} -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
        """


## Mapping: 2, spike


rule map_to_spike_by_bowtie2:
    input:
        r1="bowtie2_mapping/{rn}_contamination.1.fq",
        r2="bowtie2_mapping/{rn}_contamination.2.fq",
    output:
        sam=temp("bowtie2_mapping/{rn}_spike.sam"),
        un1=temp("bowtie2_mapping/{rn}_spike.1.fq"),
        un2=temp("bowtie2_mapping/{rn}_spike.2.fq"),
        report="bowtie2_mapping/{rn}_spike.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: os.path.join(ref_dir, REF["spike"]["bt2"]),
        un="bowtie2_mapping/{rn}_spike.fq",
    threads: 24
    resources:
        mem="96G",
    shell:
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc {params.un} -x {params.ref_bowtie2} -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
        """


## Mapping: 3, rRNA


rule map_to_rRNA_by_bowtie2:
    input:
        r1="bowtie2_mapping/{rn}_spike.1.fq",
        r2="bowtie2_mapping/{rn}_spike.2.fq",
    output:
        sam=temp("bowtie2_mapping/{rn}_rRNA.sam"),
        un1=temp("bowtie2_mapping/{rn}_rRNA.1.fq"),
        un2=temp("bowtie2_mapping/{rn}_rRNA.2.fq"),
        report="bowtie2_mapping/{rn}_rRNA.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: os.path.join(ref_dir, REF["rRNA"]["bt2"]),
        un="bowtie2_mapping/{rn}_rRNA.fq",
    threads: 24
    resources:
        mem="96G",
    shell:
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc {params.un} -x {params.ref_bowtie2} -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
        """


## Mapping: 4, smallRNA


rule map_to_smallRNA_by_bowtie2:
    input:
        r1="bowtie2_mapping/{rn}_rRNA.1.fq",
        r2="bowtie2_mapping/{rn}_rRNA.2.fq",
    output:
        sam=temp("bowtie2_mapping/{rn}_smallRNA.sam"),
        un1=temp("bowtie2_mapping/{rn}_smallRNA.1.fq"),
        un2=temp("bowtie2_mapping/{rn}_smallRNA.2.fq"),
        report="bowtie2_mapping/{rn}_smallRNA.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: os.path.join(ref_dir, REF["smallRNA"]["bt2"]),
        un="bowtie2_mapping/{rn}_smallRNA.fq",
    threads: 24
    resources:
        mem="96G",
    shell:
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc {params.un} -x {params.ref_bowtie2} -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
        """


# post mapping


rule sort_and_filter_bam_bowtie2:
    input:
        "bowtie2_mapping/{rn}_{reftype}.sam",
    output:
        "run_mapping/{rn}_{reftype}.bam",
    wildcard_constraints:
        reftype="contamination|spike|rRNA|smallRNA",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    resources:
        mem="32G",
    shell:
        """
        {params.path_samtools} sort -@ {threads} --input-fmt-option 'filter=[NM]<=10' -m 2G -O BAM -o {output} {input}
        """


## Mapping: 5, genome


rule map_to_genome_by_star:
    input:
        "bowtie2_mapping/{rn}_smallRNA.1.fq",
        "bowtie2_mapping/{rn}_smallRNA.2.fq",
    output:
        bam="run_mapping/{rn}_genome.bam",
        log="star_mapping/{rn}_genome_Log.final.out",
        fq_1=temp("star_mapping/{rn}_genome_Unmapped.out.mate1"),
        fq_2=temp("star_mapping/{rn}_genome_Unmapped.out.mate2"),
    params:
        star_ref=lambda wildcards: os.path.join(ref_dir, REF["genome"]["star"]),
        star_path=config["path"]["star"],
        output_pre="star_mapping/{rn}_genome_",
        bam="star_mapping/{rn}_genome_Aligned.sortedByCoord.out.bam",
    wildcard_constraints:
        reftype="genome",
    threads: 24
    resources:
        mem="96G",
    shell:
        """
        ulimit -n 20000
        {params.star_path} \
          --runThreadN {threads} \
          --genomeDir {params.star_ref} \
          --readFilesIn {input} \
          --alignEndsType Local \
          --outFilterMatchNminOverLread 0.66 \
          --outFilterMatchNmin 15 \
          --outFilterMismatchNmax 5 \
          --outFilterMismatchNoverLmax 0.2 \
          --outFilterMultimapNmax 50 \
          --outSAMmultNmax -1 \
          --outReadsUnmapped Fastx \
          --outSAMattrRGline ID:{wildcards.rn} SM:{wildcards.rn} LB:RNA PL:Illumina PU:SE \
          --outSAMattributes NH HI AS nM NM MD jM jI MC \
          --limitBAMsortRAM 8000000000 \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix {params.output_pre}
        mv {params.bam} {output.bam}
        """


rule compress_star_unmap:
    input:
        "star_mapping/{rn}_genome_Unmapped.out.mate{rd}",
    output:
        "star_mapping/{rn}_genome_Unmapped.out.mate{rd}.fq.gz",
    threads: 12
    resources:
        mem="48G",
    shell:
        """
        bgzip -@ {threads} -l 9 -c {input} > {output}
        """


rule report_star:
    input:
        expand(
            "star_mapping/{rn}_{reftype}_Log.final.out",
            rn=run_ids,
            reftype=["genome"],
        ),
    output:
        "quality_control/report_star.html",
    shell:
        "multiqc -f -m star -n {output} {input}"


# unmapped-QC


rule falco_unmap:
    input:
        "star_mapping/{rn}_genome_Unmapped.out.mate{rd}.fq.gz",
    output:
        html="quality_control/falco_unmap/{rn}_R{rd}/fastqc_report.html",
        text="quality_control/falco_unmap/{rn}_R{rd}/fastqc_data.txt",
        summary="quality_control/falco_unmap/{rn}_R{rd}/summary.txt",
    params:
        "quality_control/falco_unmap/{rn}_R{rd}",
    shell:
        "falco -o {params} {input}"


rule report_falco_unmap:
    input:
        reports=expand(
            "quality_control/falco_unmap/{rn}_{rd}/fastqc_data.txt",
            rn=run_ids,
            rd=read_ids,
        ),
    output:
        "quality_control/report_falco_unmap.html",
    shell:
        "multiqc -f -m fastqc -n {output} {input.reports}"


################################################################################

# combine mapping results (multi run)


rule combine_runs:
    input:
        lambda wildcards: [
            f"run_mapping/{r}_{wildcards.reftype}.bam"
            for r in sample2run[wildcards.sample]
        ],
    output:
        temp("combined_mapping/{sample}_{reftype}.bam"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    resources:
        mem="16G",
    shell:
        "{params.path_samtools} merge -@ {threads} -o {output} {input}"


## stat reads


rule count_cutadapt_reads:
    input:
        lambda wildcards: [
            f"cut_adapter/{r}_step{j}.report"
            for r in sample2run[wildcards.sample]
            for j in [1, 2]
        ],
    output:
        "stat_reads/{sample}.tsv",
    params:
        py=os.path.join(src_dir, "parse_cutadapt_report.py"),
    shell:
        """
        {params.py} {input} >{output}
        """


## stat mapping


rule separate_genome_combined:
    input:
        "combined_mapping/{sample}_genome.bam",
    output:
        uniq=temp("separate_genome_combined/{sample}_{reftype}_unique.bam"),
        multi=temp("separate_genome_combined/{sample}_{reftype}_multi.bam"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        """
        {params.path_samtools} view -@ {threads} --input-fmt-option 'filter=[NH]==1' -h {input} -O BAM -o {output.uniq}
        {params.path_samtools} view -@ {threads} --input-fmt-option 'filter=[NH]>1' -h {input} -O BAM -o {output.multi}
        """


rule stat_mapping:
    input:
        bam=lambda wildcards: [
            f"combined_mapping/{wildcards.sample}_{reftype}.bam"
            for reftype in ["contamination", "spike", "rRNA", "smallRNA"]
        ]
        + [
            f"separate_genome_combined/{wildcards.sample}_{reftype_maptype}.bam"
            for reftype_maptype in ["genome_unique", "genome_multi"]
        ],
    output:
        tsv="stat_mapping/{sample}.tsv",
    params:
        path_samtools=config["path"]["samtools"],
        ref=[
            "contamination",
            "spike",
            "rRNA",
            "smallRNA",
            "genome_unique",
            "genome_multi",
        ],
    threads: 2
    resources:
        mem="8G",
    shell:
        # {params.path_samtools} flagstats -@ {threads} -O tsv $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}$3 == "mapped"{{t=$1}}$3 == "primary mapped"{{p=$1}}END{{print ref,p; if(t > p)print ref"_multi",t-p}}' >> {output}
        # echo -e "Input\\t"$(echo $(zcat {input.fq}|wc -l)/4|bc) > {output}
        # Already divided by 2 (paired-end)
        """
        paste <(echo {params.ref} |  tr " " "\n") <(echo {input.bam} |  tr " " "\n") | while read ref file; do
            {params.path_samtools} flagstats -@ {threads} -O tsv $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}$3 == "primary mapped"{{print ref,int($1/2)}}' >> {output}
        done
        """


## drop duplicates


rule drop_duplicates:
    input:
        bam="combined_mapping/{sample}_{reftype}.bam",
    output:
        bam="drop_duplicates/{sample}_{reftype}.bam",
        log="drop_duplicates/{sample}_{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
    threads: 8
    resources:
        mem="64G",
    shell:
        """
        java -server -Xms4G -Xmx64G -Xss100M -jar {params.path_umicollapse} bam \
            --two-pass -i {input.bam} -o {output.bam}  >{output.log}
        """


rule dedup_bam_index:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
    output:
        "drop_duplicates/{sample}_{reftype}.bam.bai",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} index -@ {threads} {input}"


## stat read count after drop duplicate


rule separate_genome_dedup:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
    output:
        uniq=temp("separate_dedup/{sample}_{reftype}_unique.bam"),
        multi=temp("separate_dedup/{sample}_{reftype}_multi.bam"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        """
        {params.path_samtools} view -@ {threads} --input-fmt-option 'filter=[NH]==1' -h {input} -O BAM -o {output.uniq}
        {params.path_samtools} view -@ {threads} --input-fmt-option 'filter=[NH]>1' -h {input} -O BAM -o {output.multi}
        """


rule stat_dedup:
    input:
        bam=lambda wildcards: [
            f"drop_duplicates/{wildcards.sample}_{reftype}.bam"
            for reftype in ["contamination", "spike", "rRNA", "smallRNA"]
        ]
        + [
            f"separate_dedup/{wildcards.sample}_{reftype_maptype}.bam"
            for reftype_maptype in ["genome_unique", "genome_multi"]
        ],
    output:
        tsv="stat_dedup/{sample}.tsv",
    params:
        path_samtools=config["path"]["samtools"],
        ref=[
            "contamination",
            "spike",
            "rRNA",
            "smallRNA",
            "genome_unique",
            "genome_multi",
        ],
    threads: 2
    resources:
        mem="8G",
    shell:
        """
        paste <(echo {params.ref} |  tr " " "\n") <(echo {input.bam} |  tr " " "\n") | while read ref file; do
            {params.path_samtools} flagstats -@ {threads} -O tsv $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}$3 == "primary mapped"{{print ref,int($1/2)}}' >> {output}
        done
        """


## RNAseq QC


rule rnaseq_qc:
    input:
        "drop_duplicates/{sample}_genome.bam",
    output:
        "rnaseq_qc/{sample}.metrics.tsv",
    params:
        gtf=lambda wildcards: os.path.join(ref_dir, REF["genome"]["gtf_collapse"]),
        path_rnaseqc=config["path"]["rnaseqc"],
        outdir="rnaseq_qc",
    shell:
        """
        {params.path_rnaseqc} {params.gtf} {input} {params.outdir} -s {wildcards.sample} --coverage -v
        """


rule report_rnaseqc:
    input:
        expand(
            "rnaseq_qc/{sample}.metrics.tsv",
            sample=sample_ids,
        ),
    output:
        "quality_control/report_rnaseqc.html",
    shell:
        "multiqc -f -m rna_seqc -n {output} {input}"


################################################################################

## redo spikin alignment by blast


rule merge_map:
    input:
        "drop_duplicates/{sample}_spike.bam",
    output:
        r1=temp("spike_reads_tmp/ismap/{sample}_R1.fq"),
        r2=temp("spike_reads_tmp/ismap/{sample}_R2.fq"),
        fq=temp("spike_reads_tmp/ismap/{sample}.fq"),
    threads: 12
    resources:
        mem="36G",
    shell:
        """
        samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input} 
        fastp -w {threads} -j /dev/null -h /dev/null -i {output.r1} -I {output.r2} -m --merged_out {output.fq}
        """


rule merge_unmap:
    input:
        r1="star_mapping/{rn}_genome_Unmapped.out.mate1.fq.gz",
        r2="star_mapping/{rn}_genome_Unmapped.out.mate2.fq.gz",
    output:
        temp("spike_reads_by_run/unmap/{rn}.fq"),
    threads: 12
    resources:
        mem="36G",
    shell:
        """
        fastp -w {threads} -j /dev/null -h /dev/null -i {input.r1} -I {input.r2} -m --stdout | {{ grep --no-group-separator -A 2 -B 1 -P "CTAGAATTACACCA|TGGTGTAATTCTAG" || true; }} > {output}
        """


rule combined_unmap:
    input:
        lambda wildcards: [
            f"spike_reads_by_run/unmap/{r}.fq" for r in sample2run[wildcards.sample]
        ],
    output:
        temp("spike_reads_tmp/unmap/{sample}.fq"),
    threads: 12
    resources:
        mem="36G",
    shell:
        """
        cat {input} > {output}
        """


rule spike_fasta:
    input:
        "spike_reads_tmp/ismap/{sample}.fq",
        "spike_reads_tmp/unmap/{sample}.fq",
    output:
        fq=temp("spike_reads/{sample}.fq"),
        fa=temp("spike_reads/{sample}.fa"),
    shell:
        """
        cat {input} | tee {output.fq} | seqtk seq -a > {output.fa}
        """


rule map_to_spikin_by_blastn:
    input:
        "spike_reads/{sample}.fa",
    output:
        temp("spike_aligned/{sample}.xml"),
    params:
        ref_blast=os.path.join(ref_dir, REF["spikeN"]["blast"])
    threads: 24
    resources:
        mem="48G",
    shell:
        """
        blastn -num_threads {threads} -max_target_seqs 1 -db {params.ref_blast} -query {input} -outfmt 5 > {output}
        """


rule blastn_to_bam:
    input:
        xml="spike_aligned/{sample}.xml",
        fq="spike_reads/{sample}.fq",
    output:
        temp("spike_aligned_tmp/{sample}.unsort.bam"),
    params:
        path_blast2bam=config["path"]["blast2bam"],
        ref_fa=os.path.join(ref_dir, REF["spikeN"]["fa"])
    threads: 2
    resources:
        mem="8G",
    shell:
        """
        {params.path_blast2bam} {input.xml} {params.ref_fa} {input.fq} | \
            samtools calmd -@ {threads} --input-fmt-option 'filter=pos < 10 && pos + qlen > 33 && !flag.unmap' --output-fmt BAM - {params.ref_fa} 2>/dev/null > {output}
        """


rule blastn_bam_sort:
    input:
        "spike_aligned_tmp/{sample}.unsort.bam",
    output:
        "spike_aligned/{sample}.bam",
    threads: 4
    resources:
        mem="16G",
    shell:
        """
        samtools sort -@ {threads} -m 12G --write-index {input} -o {output}
        """


rule call_mutation_of_spike:
    input:
        "spike_aligned/{sample}.bam",
    output:
        "spike_aligned/{sample}.tsv.gz",
    params:
        py=os.path.join(src_dir, "call_spike_mutation.py"),
        header="\t".join(["ref", "motif", "base", "count"]),
    threads: 4
    resources:
        mem="12G",
    shell:
        """
        (
          echo {params.header:q}
          {params.py} {input} | awk '{{ t[$0]++ }} END{{ for (i in t) print t[i],i }}' | awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$4,$1}}'
        ) | bgzip -@ {threads} -l 9 >{output}
        """


## Count read


rule count_genome_multiple:
    input:
        expand(
            "drop_duplicates/{sample}_{reftype}.bam",
            sample=sample_ids,
            reftype=["genome"],
        ),
    output:
        "count_reads/genome_single.count",
    params:
        gtf=lambda wildcards: os.path.join(ref_dir, REF["genome"]["gtf"]),
        path_featureCounts=config["path"]["featureCounts"],
    threads: 32
    resources:
        mem="42G",
    shell:
        "{params.path_featureCounts} -T {threads} -O --largestOverlap -t exon -g gene_name -a {params.gtf} -o {output} {input}"


## call mutation
## (call mutation by sample gruop, too many samples)

## calcuate sequence coverage of all A bases on RNA


rule get_covered_positions_by_group:
    input:
        bam="drop_duplicates/{sample}_{reftype}.bam",
        bai="drop_duplicates/{sample}_{reftype}.bam.bai",
    output:
        "count_depth_by_sample/{sample}_{reftype}.tsv.gz",
    params:
        path_samtools=config["path"]["samtools"],
        ref=lambda wildcards: os.path.join(ref_dir, REF[wildcards.reftype]["fa"]),
    threads: 1
    resources:
        mem="4G",
    # m6a lib is reverse
    # do not use samtools depth, it can not output the reference base
    shell:
        """
        (
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 99 == 99 || flag & 147 == 147)' --no-output-del --no-output-ins --no-output-ends -d 0 -Q 10 -f {params.ref} {input.bam} | awk '$3 == "T" || $3 == "t"' | sed 's/\\t/\\t-\\t/3'
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 83 == 83 || flag & 163 == 163)' --no-output-del --no-output-ins --no-output-ends -d 0 -Q 10 -f {params.ref} {input.bam} | awk '$3 == "A" || $3 == "a"' | sed 's/\\t/\\t+\\t/3'
        ) | cut -f 1-5 |  bgzip -@ {threads} -l 9 >{output}
        """


## prefilter position by RNAseqmut


rule get_mutation_positions_by_group:
    input:
        bam="drop_duplicates/{sample}_{reftype}.bam",
        bai="drop_duplicates/{sample}_{reftype}.bam.bai",
    output:
        temp("filter_positions_sep/{sample}_{reftype}.bed"),
    params:
        path_bedtools=config["path"]["bedtools"],
        path_rnaseqmut=config["path"]["rnaseqmut"],
        ref=lambda wildcards: os.path.join(ref_dir, REF[wildcards.reftype]["fa"]),
        max_mismatch=5,
        # report_indel=" -d ",
        report_indel="",
    threads: 1
    resources:
        mem="4G",
    shell:
        ## rnaseqmut report format:
        ## chr1    4776457 A       C       135     56      1       0
        ## This mutation occurs at chr1:4776457, and it is A to C mutation. The next 4 numbers are 
        ## the number of reference reads (forward), reference reads (backward), alternative reads (forward) and alternative reads (backward).
        ##
        ## IMPORTANT: rnaseqmut can not deal with overlap paired end reads correctly, so, just prefilter only...
        ## IMPORTANT: rnaseqmut can not infer refbase correctly, reference fasta is required...
        ## max mismatch == 5
        ## at lease one mut event
        ## The output file in the rule has 5 fields: chr, start, end, ref, mut_num, depth_num.
        """
        {params.path_rnaseqmut} -s {params.max_mismatch} -m 0 -i 1 -k {params.report_indel} -r {params.ref} {input.bam} | awk 'BEGIN{{FS="\\t";OFS="\\t"}}($3=="A" || $3 == "T"){{print $1,$2-1,$2,$3,$7+$8,$5+$6+$7+$8}}' | {params.path_bedtools} groupby -g 1,2,3,4 -c 5,6 -o sum,sum > {output}
        """


rule combine_positions_by_group:
    input:
        lambda wildcards: [
            f"filter_positions_sep/{sample}_{wildcards.reftype}.bed"
            for sample in group2sample[wildcards.group]
            if "treat" in sample
        ],
    output:
        temp("filter_positions_by_group/{group}_{reftype}_A.site"),
        temp("filter_positions_by_group/{group}_{reftype}_T.site"),
    params:
        path_bedtools=config["path"]["bedtools"],
        out_prefix="filter_positions_by_group/{group}_{reftype}_",
    shell:
        ## The input file in the rule has 5 fields: chr, start, end, ref, mut_num, depth_num.
        ## sum mut >= 2
        ## max depth >= 5
        """
        cat {input} | sort -k1,1 -k2,2n | {params.path_bedtools} groupby -i - -g 1,2,3,4 -c 5,6 -o sum,max | awk '$5>=2 && $6 >=5 {{print $1"\\t"$2"\\t"$3 > "{params.out_prefix}"$4".site"}}'
        """


rule merge_position_into_region_by_group:
    input:
        "filter_positions_by_group/{group}_{reftype}_{refbase}.site",
    output:
        "filter_positions_by_group/{group}_{reftype}_{refbase}.bed",
    params:
        path_bedtools=config["path"]["bedtools"],
    shell:
        # bedtools merge -d 50 -i {input} | bedtools makewindows -b - -w 50 | 
        """
        {params.path_bedtools} merge -i {input} > {output}
        """


## For genome, which is too large. Split it !!!


rule prepare_splitted_bed:
    input:
        "filter_positions_by_group/{group}_{reftype}_{refbase}.bed",
    output:
        [
            temp(
                f"filter_positions_by_group_split/{{group}}_{{reftype}}_{{refbase}}.{bed_index:05}.bed"
            )
            for bed_index in range(1, 1 + bed_number)
        ],
    params:
        pre="filter_positions_by_group_split/{group}_{reftype}_{refbase}",
        path_bedtools=config["path"]["bedtools"],
    threads: 2
    resources:
        mem="8G",
    shell:
        """
        {params.path_bedtools} split -i {input} -n {bed_number} -p {params.pre} -a size
        """


rule count_bases_by_group_splited:
    input:
        bam=lambda wildcards: [
            f"drop_duplicates/{sample}_{wildcards.reftype}.bam"
            for sample in group2sample[wildcards.group]
        ],
        bai=lambda wildcards: [
            f"drop_duplicates/{sample}_{wildcards.reftype}.bam.bai"
            for sample in group2sample[wildcards.group]
        ],
        bedA="filter_positions_by_group_split/{group}_{reftype}_A.{bed_index}.bed",
        bedT="filter_positions_by_group_split/{group}_{reftype}_T.{bed_index}.bed",
    output:
        temp("pileup_bases_by_group_split/{group}_{reftype}_{bed_index}.tsv"),
    params:
        path_samtools=config["path"]["samtools"],
        path_cpup=config["path"]["cpup"],
        ref=lambda wildcards: os.path.join(ref_dir, REF[wildcards.reftype]["fa"]),
    threads: 2
    resources:
        mem="16G",
    # m6a lib is reverse
    shell:
        """
        (
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 99 == 99 || flag & 147 == 147)' -d 0 -Q 10 --reverse-del -f {params.ref} -l {input.bedT} {input.bam} 2>/dev/null | {params.path_cpup} -H -S -f mut:2 | sed 's/\\t/\\t-\\t/3'
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 83 == 83 || flag & 163 == 163)' -d 0 -Q 10 --reverse-del -f {params.ref} -l {input.bedA} {input.bam} 2>/dev/null | {params.path_cpup} -H -S -f mut:2 | sed 's/\\t/\\t+\\t/3'
        ) >{output}
        """


rule count_bases_combine:
    input:
        lambda wildcards: [
            f"pileup_bases_by_group_split/{wildcards.group}_{wildcards.reftype}_{bed_index:05}.tsv"
            for bed_index in range(1, 1 + bed_number)
        ],
    output:
        "pileup_bases_by_group/{group}_{reftype}.tsv.gz",
    params:
        header=lambda wildcards: "\t".join(
            ["chr", "pos", "ref_base", "strand"] + group2sample[wildcards.group]
        ),
    threads: 8
    resources:
        mem="42G",
    shell:
        """
        (
            echo {params.header:q}
            cat {input} | sort -S 40G -k1,1 -k2,2n
        ) | bgzip -@ {threads} -l 9 > {output}
        """


rule calculate_sites:
    input:
        "pileup_bases_by_group/{group}_{reftype}.tsv.gz",
    output:
        "pileup_calculate_by_group/{group}_{reftype}.tsv.gz",
    resources:
        mem="12G",
    params:
        py=os.path.join(src_dir, "calculate_sites.py"),
    shell:
        """
        {params.py} {input} {output}
        """


rule filter_sites:
    input:
        "pileup_calculate_by_group/{group}_{reftype}.tsv.gz",
    output:
        "pileup_filtered_by_group/{group}_{reftype}.tsv.gz",
    params:
        py=os.path.join(src_dir, "stat_and_filter_sites.py"),
        fa=lambda wildcards: os.path.join(ref_dir, REF[wildcards.reftype]["fa"]),
    resources:
        mem="12G",
    shell:
        """
        {params.py} {params.fa} {input} {output}
        """
