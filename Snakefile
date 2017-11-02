rule all:
    input:
        expand('output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R{direction}_{lane}_trim.fastq',
                direction=[1,2],
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

rule trimmomatic_lab_sample:
    input:
        'input/lab_sample/{filename}_R{direction}_00{lane}.fastq',
        adapters='input/TruSeq3-PE-2.fa'
    output:
        'output/trimmomatic_lab_sample/{filename}_R{direction}_{lane}_trim.fastq'
    shell: '''
        conda load trimmomatic
        wget https://anonscm.debian.org/cgit/debian-med/trimmomatic.git/plain/adapters/TruSeq3-PE-2.fa
        mv TruSeq3-PE-2.fa input/
        
        trimmomatic PE -phred33 \
        expand('input/lab_sample/lab_sample_39872_GTGAAA_L002_R1_00{lane}.fastq,
                lane=range(1,7)
        input/lab_sample/lab_sample_39872_GTGAAA_L002_R2_{lane}.fastq \        
        output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_R1_{lane}_trim.fastq \
        output/trimmomatic_lab_sample/s1_se \
        output/trimmomatic_lab_sample/lab_sample_39872_GTGAAA_L002_{lane}_trim.fastq) \
        output/trimmomatic_lab_sample/s2_se
        ILLUMINACLIP:{adapters}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:25
        gzip -9c output/trimmomatic_lab_sample/s{1,2}_se >> output/trimmomatic_lab_sample/orphans.fq.gz
        rm -f output/trimmomatic_lab_sample s{1,2}_se ''' 
