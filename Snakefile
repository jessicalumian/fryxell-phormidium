rule all:
    input:
        expand('output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}_paired_trim.fastq',
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
        forward='input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_{lane}.fastq',
        reverse='input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_{lane}.fastq',
        adapters='input/TruSeq3-PE-2.fa'
    output:
        forward_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_{lane}_paired_trim.fastq',
        forward_unpaired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_{lane}_unpaired_trim.fastq',
        reverse_paired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_{lane}_paired_trim.fastq',
        reverse_unpaired='output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R2_{lane}_unpaired_trim.fastq'
    message:
        'Trimming Illumina adapters from {input.forward} and {input.reverse}'
    conda:
        'envs/trimmomatic.yaml'
    shell: '''
        trimmomatic PE {input.forward} {input.reverse} {output.forward_paired} \
        {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 '''

