rule all:
    input:
        interleave_out=expand('output/interleave_lab_sample/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                               lane=range(1,7))

rule fastqc_reads:
    input:
        'input/{sample}/{filename}_R{direction}_{lane}.fastq',
    output:
        'output/fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.html',
        'output/fastqc_reads/{sample}/{filename}_R{direction}_{lane}_fastqc.zip',
    shell: '''
        module load fastqc/0.11.5
        fastqc -o `dirname {output[0]}` {input}
    '''
rule download_adapters:
    output:
        adapters='input/TruSeq3-PE-2.fa'
    shell: '''
        wget https://anonscm.debian.org/cgit/debian-med/trimmomatic.git/plain/adapters/TruSeq3-PE-2.fa
        mv TruSeq3-PE-2.fa input/ '''

rule trimmomatic_lab_sample:
    input:
        forward=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq',
                        lane=range(1,7)),
        reverse=expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}.fastq',
                        lane=range(1,7)),
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
        forward_paired=expand('output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
                               lane=range(1,7)),
        reverse_paired=expand('output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_00{lane}_paired_trim.fastq',
                               lane=range(1,7))
    output:
        interleave_out=expand('output/interleave_lab_sample/lab_sample_39872_GTGAAA_L002_00{lane}_paired_trim_interleaved.fastq',
                               lane=range(1,7))
    message:
        'Interleaving {input.forward_paired} and {input.reverse_paired}'
    conda:
        'envs/interleave.yaml'
    shell: '''
        interleave-reads.py {input.forward_paired} {input.reverse_paired} {output.interleave_out} '''
