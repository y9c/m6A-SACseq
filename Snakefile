import sys
from collections import defaultdict


workdir: config["workdir"]


src_dir = config["srcdir"]
tmp_dir = config["tmpdir"]
data_dir = os.path.dirname(workflow.configfiles[-1])

if "samples" not in config:
    sys.exit("`samples` is not defined in config file!")

if "references" not in config:
    sys.exit("`references` is not defined in config file!")

REF = config["references"]
for k, v in REF.items():
    for k2, v2 in v.items():
        v2 = os.path.expanduser(v2)
        if not os.path.isabs(v2):
            REF[k][k2] = os.path.join(data_dir, v2)
        else:
            REF[k][k2] = v2


REFTYPES = ["contamination", "genes"] if "contamination" in REF else ["genes"]


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
                run2file[s + f"-run{i}"] = {
                    x: os.path.expanduser(y) for x, y in r.items()
                }
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
        "count_reads/genome_single.count",
        expand(
            "count_depth_by_sample/{sample}_{reftype}.tsv.gz",
            sample=sample_ids,
            reftype=["genome"],
        ),
        expand(
            "pileup_filtered_by_group/{group}_{reftype}.tsv.gz",
            group=group_ids,
            reftype=["genes", "genome"],
        ),


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
        ref_bowit2=REF["contamination"]["bt2"] if "contamination" in REF else "",
        un="bowtie2_mapping/{rn}_contamination.fq",
    threads: 24
    shell:
        # The condition is not that stringent
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --no-unal --end-to-end --fast --un-conc {params.un} -x {params.ref_bowtie2} -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
        """


## Mapping: 2, spike + rRNA + tRNA + snRNA, etc.


rule map_to_genes_by_bowtie2:
    input:
        [
            "bowtie2_mapping/{rn}_contamination.1.fq",
            "bowtie2_mapping/{rn}_contamination.2.fq",
        ]
        if "contamination" in REF
        else ["cut_adapter/{rn}_cut_R1.fq.gz", "cut_adapter/{rn}_cut_R2.fq.gz"],
    output:
        sam=temp("bowtie2_mapping/{rn}_genes.sam"),
        un1=temp("bowtie2_mapping/{rn}_genes.1.fq"),
        un2=temp("bowtie2_mapping/{rn}_genes.2.fq"),
        report="bowtie2_mapping/{rn}_genes.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=REF["genes"]["bt2"],
        un="bowtie2_mapping/{rn}_genes.fq",
    threads: 24
    shell:
        """
        export LC_ALL="C"
        {params.path_bowtie2} -p {threads} --nofw --all --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc {params.un} -x {params.ref_bowtie2} -1 {input[0]} -2 {input[1]} > {output.sam} 2> >(tee {output.report} >&2)
        """


# post mapping


rule sort_and_filter_bam_bowtie2:
    input:
        "bowtie2_mapping/{rn}_{reftype}.sam",
    output:
        "run_mapping/{rn}_{reftype}.bam",
    wildcard_constraints:
        reftype="contamination|genes",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    shell:
        """
        {params.path_samtools} sort -@ {threads} --input-fmt-option 'filter=[NM]<=10' -m 2G -O BAM -o {output} {input}
        """


## Mapping: 5, genome


rule map_to_genome_by_star:
    input:
        "bowtie2_mapping/{rn}_genes.1.fq",
        "bowtie2_mapping/{rn}_genes.2.fq",
    output:
        bam="run_mapping/{rn}_genome.bam",
        log="star_mapping/{rn}_genome_Log.final.out",
        fq_1=temp("star_mapping/{rn}_genome_Unmapped.out.mate1"),
        fq_2=temp("star_mapping/{rn}_genome_Unmapped.out.mate2"),
    params:
        star_ref=REF["genome"]["star"],
        star_path=config["path"]["star"],
        output_pre="star_mapping/{rn}_genome_",
        bam="star_mapping/{rn}_genome_Aligned.sortedByCoord.out.bam",
    wildcard_constraints:
        reftype="genome",
    threads: 24
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
            for reftype in REFTYPES
        ]
        + [
            f"separate_genome_combined/{wildcards.sample}_{reftype_maptype}.bam"
            for reftype_maptype in ["genome_unique", "genome_multi"]
        ],
    output:
        tsv="stat_mapping/{sample}.tsv",
    params:
        path_samtools=config["path"]["samtools"],
        ref=REFTYPES
        + [
            "genome_unique",
            "genome_multi",
        ],
    threads: 2
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
    shell:
        """
        export TMPDIR={tmp_dir}
        java -server -Xms4G -Xmx64G -Xss100M -Djava.io.tmpdir={tmp_dir} -jar {params.path_umicollapse} bam \
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
            f"drop_duplicates/{wildcards.sample}_{reftype}.bam" for reftype in REFTYPES
        ]
        + [
            f"separate_dedup/{wildcards.sample}_{reftype_maptype}.bam"
            for reftype_maptype in ["genome_unique", "genome_multi"]
        ],
    output:
        tsv="stat_dedup/{sample}.tsv",
    params:
        path_samtools=config["path"]["samtools"],
        ref=REFTYPES
        + [
            "genome_unique",
            "genome_multi",
        ],
    threads: 2
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
        gtf=REF["genome"]["gtf"],
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
        gtf=REF["genome"]["gtf"],
        path_featureCounts=config["path"]["featureCounts"],
    threads: 32
    shell:
        "{params.path_featureCounts} -T {threads} -O --largestOverlap -t exon -g gene_name -a {params.gtf} -o {output} {input}"


## call mutation
## (call mutation by sample gruop, too many samples)


## calcuate sequence coverage of all A bases on RNA


rule get_covered_positions_by_group:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
    output:
        "count_depth_by_sample/{sample}_{reftype}.tsv.gz",
    params:
        path_samtools=config["path"]["samtools"],
        ref=lambda wildcards: REF[wildcards.reftype]["fa"],
    threads: 1
    # m6a lib is reverse
    # do not use samtools depth, it can not output the reference base
    shell:
        """
        (
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 99 == 99 || flag & 147 == 147)' --no-output-del --no-output-ins --no-output-ends -d 0 -Q 10 -f {params.ref} {input} | awk '$3 == "T" || $3 == "t"' | sed 's/\\t/\\t-\\t/3'
          {params.path_samtools} mpileup --input-fmt-option 'filter=(flag & 83 == 83 || flag & 163 == 163)' --no-output-del --no-output-ins --no-output-ends -d 0 -Q 10 -f {params.ref} {input} | awk '$3 == "A" || $3 == "a"' | sed 's/\\t/\\t+\\t/3'
        ) | cut -f 1-5 |  bgzip -@ {threads} -l 9 >{output}
        """


############

## prefilter position by merge all treated samples


rule merge_mutated_treated_bam:
    input:
        lambda wildcards: [
            f"drop_duplicates/{sample}_{wildcards.reftype}.bam"
            for sample in group2sample[wildcards.group]
            if "treat" in sample
        ],
    output:
        bam=temp("merged_mutated_reads_by_group/{group}_{reftype}.bam"),
        bai=temp("merged_mutated_reads_by_group/{group}_{reftype}.bam.bai"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 16
    shell:
        "{params.path_samtools} merge --write-index -@ {threads} --input-fmt-option 'filter=[NM]>0' -o {output.bam}##idx##{output.bai} {input}"


rule prefilter_positions_by_group:
    input:
        bam="merged_mutated_reads_by_group/{group}_{reftype}.bam",
        bai="merged_mutated_reads_by_group/{group}_{reftype}.bam.bai",
    output:
        "filter_positions_by_group/{group}_{reftype}_{refbase}.bed",
    params:
        ref=lambda wildcards: REF[wildcards.reftype]["fa"],
        flag=lambda wildcards: "83 163" if wildcards.refbase == "A" else "99 147",
        strand=lambda wildcards: "+" if wildcards.refbase == "A" else "-",
        path_caller=config["path"]["caller"],
    threads: 16
    shell:
        """
        {params.path_caller} -t {threads}  -i {input.bam} -r {params.ref} -b {wildcards.refbase} -f {params.flag} -d 3 -m 1 -F 3584 | \
            awk 'BEGIN{{OFS="\\t"}}{{print $1,$2-1,$2,$4"/"$5,$4/$5,"{params.strand}"}}' >{output}
        """


## Count with in prefilter positions


rule count_site_by_sample:
    input:
        bam="drop_duplicates/{sample}_{reftype}.bam",
        bai="drop_duplicates/{sample}_{reftype}.bam.bai",
        bed="filter_positions_by_group/{group}_{reftype}_{refbase}.bed",
    output:
        temp("pileup_bases_by_sample/{group}_{sample}_{reftype}_{refbase}.tsv"),
    params:
        ref=lambda wildcards: REF[wildcards.reftype]["fa"],
        flag=lambda wildcards: "83 163" if wildcards.refbase == "A" else "99 147",
        path_caller=config["path"]["caller"],
    threads: 16
    shell:
        """
        {params.path_caller} -t {threads} -i {input.bam} -s {input.bed} -r {params.ref} -b {wildcards.refbase} -f {params.flag} -d 0 -m 0 -F 3584 >{output}
        """


rule join_sites_by_group:
    input:
        lambda wildcards: expand(
            "pileup_bases_by_sample/{{group}}_{sample}_{{reftype}}_{refbase}.tsv",
            refbase=["A", "T"],
            sample=group2sample[wildcards.group],
        ),
    output:
        "pileup_calculate_by_group/{group}_{reftype}.tsv.gz",
    params:
        sample=lambda wildcards: group2sample[wildcards.group],
        py=os.path.join(src_dir, "join_samples_sites.py"),
    shell:
        """
        {params.py} -f {input} -n {params.sample} {params.sample} -o {output}
        """


rule filter_sites:
    input:
        "pileup_calculate_by_group/{group}_{reftype}.tsv.gz",
    output:
        "pileup_filtered_by_group/{group}_{reftype}.tsv.gz",
    params:
        py=os.path.join(src_dir, "filter_group_sites.py"),
        fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    shell:
        """
        {params.py} {params.fa} {input} {output}
        """
