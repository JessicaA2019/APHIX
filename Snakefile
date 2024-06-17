########################
### FIXED PARAMETERS ###
########################
reference_fasta = config.get("reference_fasta")
if not reference_fasta:
    raise RuntimeError("No reference FASTA found. Please specify 'reference_fasta' in config file")

reference_gff = config.get("reference_gff_name")

input_folder = config.get('input_fastq')
if not input_folder:
    raise RuntimeError("No input FASTQ files found. Please specify 'input_fastq' in config file")

sample_name = config.get("sample_name")

threads = config.get("threads")

if config.get("input_is_fastq"):
    polish_fa_param = ""
elif not config.get("input_is_fastq"):
    polish_fa_param = "-f"

#########################
## Optional parameters ##
#########################
minimap2_param = "-ax splice"
cluster_param = "-c "+str(config.get("cluster_size"))+" -d 5 -p 0"
polish_param = "-c "+str(config.get("cluster_size"))
gff_compare_param = "-V -M -T"
isoform_checker_param = "-n True"

########################
######### RULES ########
########################

rule all:
    input:
        #expand("{name}_isoform_analysis/mapped_reads/{name}.bam", name=sample_name),
        #expand("{name}_isoform_analysis/GFF/{name}.gff", name=sample_name),
        #expand("{name}_isoform_analysis/clustering/{name}_clusters.tsv", name=sample_name),
        #expand("{name}_isoform_analysis/clustering/{name}_cluster.gff", name=sample_name),
        #expand("{name}_isoform_analysis/polishing/{name}_polished.fa", name=sample_name),
        #expand("{name}_isoform_analysis/mapped_reads/{name}_polished.bam", name=sample_name),
        #expand("{name}_isoform_analysis/GFF/{name}_polished.gff", name=sample_name),
        #expand("{name}_isoform_analysis/gff_compare/{name}_novel_junction.tab", name=sample_name),
        expand("{name}_isoform_analysis/8_corrected_isoforms/{name}.csv", name=sample_name)

rule ref_no_header:
    input:
        REF = reference_fasta
    output:
        REF = expand("{name}_isoform_analysis/1_mapped_reads/{name}_ref_no_header.fa", name=sample_name)
    shell:
        """
        seqkit seq -i -w 0 {input.REF} > {output.REF}
        """

rule map_fq:
    input:
        FQ = input_folder,
        REF = "{name}_isoform_analysis/1_mapped_reads/{name}_ref_no_header.fa"
    params:
        minimap2_param = minimap2_param
    output:
        BAM = "{name}_isoform_analysis/1_mapped_reads/{name}.bam",
        BAI = "{name}_isoform_analysis/1_mapped_reads/{name}.bam.bai"
    threads: threads
    shell:
        """
        minimap2 {params.minimap2_param} -t {threads} --secondary=no {input.REF} {input.FQ} | samtools view -bS -F4 -h -@ {threads} - | samtools sort -@ {threads} - > {output.BAM} && samtools index -@ {threads} {output.BAM}
        """
        
rule bam2gff:
    input:
        BAM = "{name}_isoform_analysis/1_mapped_reads/{name}.bam"
    output:
        GFF = "{name}_isoform_analysis/2_GFF/{name}.gff"
    ## DO NOT ADD THREADS. It results in incorrect data.
    shell:
        """
        spliced_bam2gff -M {input.BAM} > {output.GFF}
        """

rule cluster:
    input:
        GFF = "{name}_isoform_analysis/2_GFF/{name}.gff"
    params:
        cluster_param = cluster_param
    output:
        TSV = "{name}_isoform_analysis/3_clustering/{name}_clusters.tsv",
        GFF = "{name}_isoform_analysis/3_clustering/{name}_cluster.gff"
    threads: threads
    shell:
        """
        cluster_gff {params.cluster_param} -t {threads} -a {output.TSV} {input.GFF} > {output.GFF}
        """
        
rule polish:
    input:
        TSV = "{name}_isoform_analysis/3_clustering/{name}_clusters.tsv",
        BAM = "{name}_isoform_analysis/1_mapped_reads/{name}.bam"
    params:
        polish_param = polish_param,
	polish_fa_param = polish_fa_param
    output:
        FA = "{name}_isoform_analysis/4_polishing/{name}_polished.fa"
    threads: threads
    shell:
        """
        polish_clusters -a {input.TSV} {params.polish_param} {params.polish_fa_param} -t {threads} -o {output.FA} {input.BAM}
        """
        
rule map_clusters:
    input:
        FA = "{name}_isoform_analysis/4_polishing/{name}_polished.fa",
        REF = "{name}_isoform_analysis/1_mapped_reads/{name}_ref_no_header.fa"
    params:
        minimap2_param = minimap2_param
    output:
        BAM = "{name}_isoform_analysis/5_polish_map/{name}_polished.bam",
        BAI = "{name}_isoform_analysis/5_polish_map/{name}_polished.bam.bai"
    threads: threads
    shell:
        """
        minimap2 {params.minimap2_param} -t {threads} --secondary=no {input.REF} {input.FA} | samtools view -bS -F4 -h -@ {threads} - | samtools sort -@ {threads} - > {output.BAM} && samtools index -@ {threads} {output.BAM}
        """
        
rule bam2gff_clusters:
    input:
        BAM = "{name}_isoform_analysis/5_polish_map/{name}_polished.bam"
    output:
        GFF = "{name}_isoform_analysis/6_polish_gff/{name}_polished.gff"
    ## DO NOT ADD THREADS. It results in incorrect data.
    shell:
        """
        spliced_bam2gff  -M {input.BAM} > {output.GFF}
        """
rule create_gff_ref:
    input:
        POL = expand("{name}_isoform_analysis/6_polish_gff/{name}_polished.gff", name=sample_name),
	REF = expand("{name}_isoform_analysis/1_mapped_reads/{name}_ref_no_header.fa", name=sample_name),
	GFF = "auxillary_scripts/auxillary.gff"
    output:
        GFF = reference_gff
    shell:
        """
        python3 auxillary_scripts/CLI_gff_converter.py {input.GFF} {output.GFF} {input.REF} 
        """  

rule mv_gff_ref:
    input:
        GFF = reference_gff
    output:
        GFF = expand("{name}_isoform_analysis/7_gff_compare/{ref_name}", name=sample_name, ref_name=reference_gff)
    shell:
        """
        mv {input.GFF} {output.GFF}
        """

rule gff_compare:
    input:
        GFF = expand("{name}_isoform_analysis/6_polish_gff/{name}_polished.gff", name=sample_name),
        REF = expand("{name}_isoform_analysis/7_gff_compare/{ref_name}", name=sample_name, ref_name=reference_gff)
    params:
        gff_compare_param = gff_compare_param
    output:
        TAB = "gffcmp_novel_junction.tab",
        STATS = "gffcmp.stats",
        TRACK = "gffcmp.tracking",
        LOCI = "gffcmp.loci",
        GTF = "gffcmp.annotated.gtf"
    shell:
        """
        gffcompare {params.gff_compare_param} -r {input.REF} {input.GFF} -j {output.TAB}
        """
        
rule rename_gff_compare_files:
    input:
        STATS = "gffcmp.stats",
        TRACK = "gffcmp.tracking",
        LOCI = "gffcmp.loci",
        GTF = "gffcmp.annotated.gtf",
    	TAB = "gffcmp_novel_junction.tab"
    output:
        STATS = expand("{name}_isoform_analysis/7_gff_compare/{name}.stats", name=sample_name),
        TRACK = expand("{name}_isoform_analysis/7_gff_compare/{name}.tracking", name=sample_name),
        LOCI = expand("{name}_isoform_analysis/7_gff_compare/{name}.loci", name=sample_name),
        GTF = expand("{name}_isoform_analysis/7_gff_compare/{name}.annotated.gtf", name=sample_name),
        TAB = expand("{name}_isoform_analysis/7_gff_compare/{name}_novel_junction.tab", name=sample_name)
    shell:
        """
        mv {input.STATS} {output.STATS} && mv {input.TRACK} {output.TRACK} && mv {input.LOCI} {output.LOCI} && mv {input.GTF} {output.GTF} && mv {input.TAB} {output.TAB}
        """

rule HIV_Isoform_Checker:
    input:
        GFF = expand("{name}_isoform_analysis/7_gff_compare/{name}.annotated.gtf", name=sample_name),
        REF = expand("{name}_isoform_analysis/1_mapped_reads/{name}_ref_no_header.fa", name=sample_name)
    params:
        isoform_checker_param = isoform_checker_param
    output:
    	ALT = "checked_isoforms_altered.txt",
        CSV = "checked_isoforms.csv",
        FAIL = "checked_isoforms_fail.txt",
        PASS = "checked_isoforms_pass.txt",
        ISO = "checked_isoforms_isoform_counts.csv",
        LOG = "checked_isoforms.log",
        REF = "checked_isoforms_ref_coordinates.txt",
        SS = "checked_isoforms_splice_site_usage.csv"
    shell:
        """
        HIV_Isoform_Checker {params.isoform_checker_param} {input.GFF} checked_isoforms {input.REF}
        """

rule rename_checker_files:
    input:
        ALT = "checked_isoforms_altered.txt",
        CSV = "checked_isoforms.csv",
        FAIL = "checked_isoforms_fail.txt",
        PASS = "checked_isoforms_pass.txt",
        ISO = "checked_isoforms_isoform_counts.csv",
        LOG = "checked_isoforms.log",
        REF = "checked_isoforms_ref_coordinates.txt",
        SS = "checked_isoforms_splice_site_usage.csv"
    output:
        ALT = "{name}_isoform_analysis/8_corrected_isoforms/{name}_altered.txt",
        CSV = "{name}_isoform_analysis/8_corrected_isoforms/{name}.csv",
        FAIL = "{name}_isoform_analysis/8_corrected_isoforms/{name}_fail.txt",
        PASS = "{name}_isoform_analysis/8_corrected_isoforms/{name}_pass.txt",
        ISO = "{name}_isoform_analysis/8_corrected_isoforms/{name}_isoform_counts.csv",
        LOG = "{name}_isoform_analysis/8_corrected_isoforms/{name}.log",
        REF = "{name}_isoform_analysis/8_corrected_isoforms/{name}_ref_coordinates.txt",
        SS = "{name}_isoform_analysis/8_corrected_isoforms/{name}_splice_site_usage.csv"
    shell:
        """
        mv {input.ALT} {output.ALT} && mv {input.CSV} {output.CSV} && mv {input.FAIL} {output.FAIL} && mv {input.PASS} {output.PASS} && mv {input.ISO} {output.ISO} && mv {input.LOG} {output.LOG} && mv {input.REF} {output.REF} && mv {input.SS} {output.SS}
        """
