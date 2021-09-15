##To plot profile of TF from DAP-seq

#! /bin/bash

cut -f3 ../data/highpccTF/$1 |while read fq
do
    fas=`echo $fq | tr '_' '\n'|head -n1`
    bed=`echo $1 |tr 'R' '\n'|head -n1`
    # download the fastq file
    fastq-dump --split-3 --gzip $fas

    if test `ls ${fas}*.fastq.gz|wc -l` = 2; then

        trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired ${fas}_1.fastq.gz ${fas}_2.fastq.gz
        bowtie2 -p 10 -x ../data/bowtie2index/TAIR10 -1 ${fas}_1_trimmed.fq.gz -2 ${fas}_2_trimmed.fq.gz  | samtools sort \
            -O bam -@ 11 -o - > ${fas}.sorted.bam
    else

        trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 $fas.fastq.gz
		bowtie2 -p 10 -x ../data/bowtie2index/TAIR10 -U ${fas}_trimmed.fq.gz | samtools sort -O bam -@ 11 -o - \
            > ${fas}.sorted.bam
    fi

    # make index
    samtools index ${fas}.sorted.bam

    bamCoverage -b ${fas}.sorted.bam -o ${fas}.bw
    computeMatrix scale-regions -b 1000 -a 1000 -R ../data/geneClustersBED/$bed.bed -S ${fas}.bw \
        --skipZeros -o ${fas}.gz -p 6
    plotProfile -m ${fas}.gz -out ${bed}.${fas}.png

    #remove all fastq files
    rm ${fas}*

done

