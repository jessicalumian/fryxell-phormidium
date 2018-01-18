# add quast rule (and add to rule all?)

rule all:
    input:
        'output/anvio_bowtie_build_mat/contigs_fixed/anvio-contigs.db.1'

rule fastqc_reads:
    input:
        'input/{sample}/{filename}_R{direction}_{lane}.fastq',
    output:
        'output/fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.html',
        'output/fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.zip',
    shell: '''
        module load fastqc/0.11.5
        fastqc -o `dirname {output[0]}` {input} '''

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
        forward_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
        forward_unpaired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}_unpaired_trim.fastq',
        reverse_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}_paired_trim.fastq',
        reverse_unpaired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}_unpaired_trim.fastq'
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
        forward_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
        reverse_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}_paired_trim.fastq'
    output:
        interleave_out='output/interleave_lab_sample/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq'
    message:
        'Interleaving {input.forward_paired} and {input.reverse_paired}'
    conda:
        'envs/interleave.yaml'
    shell: '''
        interleave-reads.py {input.forward_paired} {input.reverse_paired} -o {output.interleave_out} '''

rule megahit_lab:
    input: 
        expand('output/interleave_lab_sample/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                lane=range(1,7))
    output:
        'output/megahit_lab/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    params:
        input_list=lambda w, input: ','.join(input)
    shell: '''
        rm -rf output/megahit_lab
        megahit --12 {params.input_list} -o output/megahit_lab '''

rule megahit_mat:
    input:
        'input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        'output/megahit_mat/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    shell: '''
        rm -fr output/megahit_mat
        megahit --12 {input} -o output/megahit_mat '''

rule megahit_coassembly:
    input:
        expand('output/interleave_lab_sample/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                lane=range(1,7)),
        'input/mat_sample/mat_sample_104_ABC_L00_R12_0.fastq'
    output:
        'output/megahit_coassembly/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    params:
        input_list=lambda w, input: ','.join(input)
    shell: '''
        rm -fr output/megahit_coassembly
        megahit --12 {params.input_list} -o output/megahit_coassembly '''

rule anvio_reform_fasta:
    input:
        'output/megahit_{sample_type}/final.contigs.fa'
    output:
        fixed_contigs='output/anvio_reform_fasta/{sample_type}/contigs_fixed.fa',
        report='output/anvio_reform_fasta/{sample_type}/name_conversions.txt'
    conda:
        'envs/anvio.yaml'
    shell: '''
        anvi-script-reformat-fasta {input} -o {output.fixed_contigs} --min-len 2000 --simplify-names --report {output.report} '''

rule anvio_bowtie_build_mat:
    input:
        'output/anvio_reform_fasta/mat/contigs_fixed.fa'
    output:
        dynamic('output/anvio_bowtie_build_mat/contigs_fixed/anvio-contigs.db.{mat_version}')
    conda:
        'envs/anvio.yaml'
    shell:
        '''
            bowtie2-build {input} {output} '''

