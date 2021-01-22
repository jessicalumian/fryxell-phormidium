# add quast rule (and add to rule all?)

rule all:
    input:
        'output/04_anvio/mat/mat.bam-sorted.bam.bai',
        'output/04_anvio/lab/lab.bam-sorted.bam.bai',
        'output/04_anvio/coassembly/coassembly.bam-sorted.bam.bai',
        'output/04_anvio/mat/anvio-contigs.1.bt2',
        'output/04_anvio/lab/anvio-contigs.1.bt2',
        'output/04_anvio/coassembly/anvio-contigs.1.bt2',
        'output/04_anvio/mat/mat_HMMs/Rinke_et_al/genes.txt',
        'output/04_anvio/mat/mat_HMMs/Campbell_et_al/genes.txt',
        'output/04_anvio/mat/mat_HMMs/Ribosomal_RNAs/genes.txt',
        'output/04_anvio/lab/lab_HMMs/Rinke_et_al/genes.txt',
        'output/04_anvio/lab/lab_HMMs/Campbell_et_al/genes.txt',
        'output/04_anvio/lab/lab_HMMs/Ribosomal_RNAs/genes.txt',
        'output/04_anvio/coassembly/coassembly_HMMs/Rinke_et_al/genes.txt',
        'output/04_anvio/coassembly/coassembly_HMMs/Campbell_et_al/genes.txt',
        'output/04_anvio/coassembly/coassembly_HMMs/Ribosomal_RNAs/genes.txt',
        'output/04_anvio/lab/lab_ANVIO-PROFILE/PROFILE.db',
        'output/04_anvio/mat/mat_ANVIO-PROFILE/PROFILE.db',
        'output/04_anvio/coassembly/coassembly_ANVIO-PROFILE/PROFILE.db',
        'output/04_anvio/merged_mat_lab_coassembly/lab.sam',
        'output/04_anvio/merged_mat_lab_coassembly/mat.sam',
        'output/04_anvio/merged_mat_lab_coassembly/lab.bam-sorted.bam.bai',
        'output/04_anvio/merged_mat_lab_coassembly/mat.bam-sorted.bam.bai',
        'output/04_anvio/merged_mat_lab_coassembly/lab_ANVIO-PROFILE/PROFILE.db',
        'output/04_anvio/merged_mat_lab_coassembly/mat_ANVIO-PROFILE/PROFILE.db',
        'output/04_anvio/merged_mat_lab_coassembly/SAMPLES-MERGED/PROFILE.db',
        'output/04_anvio/merged_mat_lab_coassembly/SAMPLES-SUMMARY/bins_summary.txt'

rule fastqc_reads:
    input:
        'input/{sample}/{filename}_R{direction}_{lane}.fastq',
    output:
        'output/00_fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.html',
        'output/00_fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.zip',
    shell: '''
        module load fastqc/0.11.5
        fastqc -o `dirname {output[0]}` {input} '''
# adding something here
rule download_adapters:
    output:
        adapters='input/TruSeq3-PE-2.fa'
    shell: '''
        wget https://anonscm.debian.org/cgit/debian-med/trimmomatic.git/plain/adapters/TruSeq3-PE-2.fa
        mv TruSeq3-PE-2.fa input/ '''

rule trimmomatic_lab_sample:
    input:
        forward='input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq',
        reverse='input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}.fastq',
        adapters='input/TruSeq3-PE-2.fa'
    output:
        forward_paired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
        forward_unpaired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R1_00{lane}_unpaired_trim.fastq',
        reverse_paired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R2_00{lane}_paired_trim.fastq',
        reverse_unpaired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R2_00{lane}_unpaired_trim.fastq'
    message:
        'Trimming Illumina adapters from {input.forward} and {input.reverse}'
    conda:
        'envs/trimmomatic.yaml'
    shell: '''
        trimmomatic PE {input.forward} {input.reverse} {output.forward_paired} \
        {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 '''

rule interleave_lab_sample:
    input:
        forward_paired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
        reverse_paired='output/01_trimmomatic_lab/lab_sample_39872_GTGAAA_L002_R2_00{lane}_paired_trim.fastq'
    output:
        interleave_out='output/02_interleave_lab/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq'
    message:
        'Interleaving {input.forward_paired} and {input.reverse_paired}'
    conda:
        'envs/interleave.yaml'
    shell: '''
        interleave-reads.py {input.forward_paired} {input.reverse_paired} -o {output.interleave_out} '''

rule megahit_lab:
    input: 
        expand('output/02_interleave_lab/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                lane=range(1,7))
    output:
        'output/03_megahit_lab/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    params:
        input_list=lambda w, input: ','.join(input)
    shell: '''
        rm -rf output/03_megahit_lab
        megahit --12 {params.input_list} -o output/03_megahit_lab '''

rule megahit_mat:
    input:
        'input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        'output/03_megahit_mat/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    shell: '''
        rm -fr output/03_megahit_mat
        megahit --12 {input} -o output/03_megahit_mat '''

rule megahit_coassembly:
    input:
        expand('output/02_interleave_lab/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                lane=range(1,7)),
        'input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        'output/03_megahit_coassembly/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    params:
        input_list=lambda w, input: ','.join(input)
    shell: '''
        rm -fr output/03_megahit_coassembly
        megahit --12 {params.input_list} -o output/03_megahit_coassembly '''

rule anvio_reform_fasta:
    input:
        'output/03_megahit_{sample}/final.contigs.fa'
    output:
        fixed_contigs='output/04_anvio/{sample}/contigs_fixed.fa',
        report='output/04_anvio/{sample}/name_conversions.txt'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-script-reformat-fasta {input} -o {output.fixed_contigs} --min-len 2000 --simplify-names --report {output.report} '''

rule anvio_bowtie_build:
    input:
        'output/04_anvio/{sample}/contigs_fixed.fa'
    output:
        'output/04_anvio/{sample}/anvio-contigs.1.bt2'
    params:
       bt2_base='output/04_anvio/{sample}/anvio-contigs' 
    conda:
        'envs/anvio.yaml'
    shell:
        '''
            bowtie2-build {input} {params.bt2_base} '''

rule bowtie2_samtools_map_mat:
    input:
        raw_reads='input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        sam='output/04_anvio/mat/mat.sam',
        bam='output/04_anvio/mat/mat.bam'
    conda:
        'envs/anvio.yaml'
    params:
        bt2_base='output/04_anvio/mat/anvio-contigs'
    shell: '''
        bowtie2 --threads 8 -x {params.bt2_base} -U {input.raw_reads} -S {output.sam}
        samtools view -U 4 -bS {output.sam} > {output.bam} '''

rule MERGE_bowtie2_samtools_map_mat:
    input:
        raw_reads='input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        sam='output/04_anvio/merged_mat_lab_coassembly/mat.sam',
        bam='output/04_anvio/merged_mat_lab_coassembly/mat.bam'
    conda:
        'envs/anvio.yaml'
    params:
        bt2_base='output/04_anvio/coassembly/anvio-contigs'
    shell: '''
        bowtie2 --threads 8 -x {params.bt2_base} -U {input.raw_reads} -S {output.sam}
        samtools view -U 4 -bS {output.sam} > {output.bam} '''

rule MERGE_bowtie2_samtools_map_lab:
    input:
        forward=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq',
                        lane=range(1,7)),
        reverse=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}.fastq',
                        lane=range(1,7))
    output:
        sam='output/04_anvio/merged_mat_lab_coassembly/lab.sam',
        bam='output/04_anvio/merged_mat_lab_coassembly/lab.bam'
    conda:
        'envs/anvio.yaml'
    params:
        input_list=lambda w, input: ','.join(input),
        bt2_base='output/04_anvio/coassembly/anvio-contigs'
    shell: '''
        bowtie2 --threads 8 -x {params.bt2_base} -U {params.input_list} -S {output.sam}
        samtools view -U 4 -bS {output.sam} > {output.bam} '''

rule bowtie2_samtools_map_lab:
    input:
        forward=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq',
                        lane=range(1,7)),
        reverse=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}.fastq',
                        lane=range(1,7))
    output:
        sam='output/04_anvio/lab/lab.sam',
        bam='output/04_anvio/lab/lab.bam'
    conda:
        'envs/anvio.yaml'
    params:
        input_list=lambda w, input: ','.join(input),
        bt2_base='output/04_anvio/lab/anvio-contigs'
    shell: ''' 
        bowtie2 --threads 8 -x {params.bt2_base} -U {params.input_list} -S {output.sam}
        samtools view -U 4 -bS {output.sam} > {output.bam} '''

rule bowtie2_samtools_map_coassembly:
    input:
        forward=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq',
                        lane=range(1,7)),
        reverse=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}.fastq',
                        lane=range(1,7)),
        mat='input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        sam='output/04_anvio/coassembly/coassembly.sam',
        bam='output/04_anvio/coassembly/coassembly.bam'
    conda:
        'envs/anvio.yaml'
    params:
        input_list=lambda w, input: ','.join(input),
        bt2_base='output/04_anvio/coassembly/anvio-contigs'
    shell: '''
        bowtie2 --threads 8 -x {params.bt2_base} -U {params.input_list} -S {output.sam}
        samtools view -U 4 -bS {output.sam} > {output.bam} '''

rule convert_bam_anvio:
    input:
        'output/04_anvio/{sample}/{sample}.bam'
    output:
        'output/04_anvio/{sample}/{sample}.bam-sorted.bam.bai'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-init-bam {input} '''

rule MERGE_convert_bam_anvio:
    input:
        'output/04_anvio/merged_mat_lab_coassembly/{sample}.bam'
    output:
        'output/04_anvio/merged_mat_lab_coassembly/{sample}.bam-sorted.bam.bai'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-init-bam {input} '''

rule anvi_gen_contigs_database:
    input:
        'output/04_anvio/{sample}/contigs_fixed.fa'
    output:
        'output/04_anvio/{sample}/anvio_contigs.db'
    params:
        '{sample}'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-gen-contigs-database -f {input} -n {params} -o {output} '''

rule run_hmms:
    input:
        'output/04_anvio/{sample}/anvio_contigs.db'
    output:
        'output/04_anvio/{sample}/{sample}_HMMs/Rinke_et_al/genes.txt',
        'output/04_anvio/{sample}/{sample}_HMMs/Campbell_et_al/genes.txt',
        'output/04_anvio/{sample}/{sample}_HMMs/Ribosomal_RNAs/genes.txt'
    params:
        org_dir='.snakemake/conda/da60b80a/lib/python3.5/site-packages/anvio/data/hmm',
        cor_output='output/04_anvio/{sample}/{sample}_HMMs'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-run-hmms -c {input} --num-threads 28 
        mkdir -p {params.cor_output}
        cp -r {params.org_dir}/Rinke_et_al {params.cor_output} 
        cp -r {params.org_dir}/Campbell_et_al {params.cor_output}
        cp -r {params.org_dir}/Ribosomal_RNAs {params.cor_output} '''

rule anvi_profile:
    input:
        bam='output/04_anvio/{sample}/{sample}.bam-sorted.bam',
        db='output/04_anvio/{sample}/anvio_contigs.db'
    output:
        'output/04_anvio/{sample}/{sample}_ANVIO-PROFILE/AUXILIARY-DATA.h5',
        'output/04_anvio/{sample}/{sample}_ANVIO-PROFILE/PROFILE.db'
    params:
        out_dir='output/04_anvio/{sample}/{sample}_ANVIO-PROFILE',
        sam_name='{sample}'
    conda:
        'envs/anvio.yaml'
    shell: '''
        rm -fr {params.out_dir}
        anvi-profile -i {input.bam} -c {input.db} -T 28 -o {params.out_dir} -S {params.sam_name} '''

rule MERGE_anvi_profile:
    input:
        bam='output/04_anvio/merged_mat_lab_coassembly/{sample}.bam-sorted.bam',
        db='output/04_anvio/coassembly/anvio_contigs.db'
    output:
        'output/04_anvio/merged_mat_lab_coassembly/{sample}_ANVIO-PROFILE/AUXILIARY-DATA.h5',
        'output/04_anvio/merged_mat_lab_coassembly/{sample}_ANVIO-PROFILE/PROFILE.db'
    params:
        out_dir='output/04_anvio/merged_mat_lab_coassembly/{sample}_ANVIO-PROFILE',
        sam_name='{sample}'
    conda:
        'envs/anvio.yaml'
    shell: '''
        rm -fr {params.out_dir}
        anvi-profile -i {input.bam} -c {input.db} -T 28 -o {params.out_dir} -S {params.sam_name} '''

rule MERGE_lab_mat_coassembly:
    input:
        lab_prof='output/04_anvio/merged_mat_lab_coassembly/lab_ANVIO-PROFILE/PROFILE.db',
        mat_prof='output/04_anvio/merged_mat_lab_coassembly/mat_ANVIO-PROFILE/PROFILE.db',
        coassembly_db='output/04_anvio/coassembly/anvio_contigs.db'
    output:
        'output/04_anvio/merged_mat_lab_coassembly/SAMPLES-MERGED/PROFILE.db'
    params:
        output_dir='output/04_anvio/merged_mat_lab_coassembly/SAMPLES-MERGED'
    conda:
        'envs/anvio.yaml'
    shell: '''
        rm -fr {params.output_dir}    
        anvi-merge {input.lab_prof} {input.mat_prof} -o {params.output_dir} -c {input.coassembly_db} --enforce-hierarchical-clustering '''

rule MERGE_anvi_summarize:
    input:
        merge_prof='output/04_anvio/merged_mat_lab_coassembly/SAMPLES-MERGED/PROFILE.db',
        coassembly_db='output/04_anvio/coassembly/anvio_contigs.db'
    output:
        'output/04_anvio/merged_mat_lab_coassembly/SAMPLES-SUMMARY/general_bins_summary.txt'
    params:
        output_dir='output/04_anvio/merged_mat_lab_coassembly/SAMPLES-SUMMARY'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-summarize -p {input.merge_prof} -c {input.coassembly_db} -o {params.output_dir} -C CONCOCT '''
