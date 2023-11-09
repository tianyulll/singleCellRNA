#!/bin/bash

'''
tianyu.lu@pennmedicine.upenn.edu

Running Cellranger in parallel on AWS.
Creates the merged fastq files from the chunk.
Rename to meet Cellranger input requirements.
Remove the merged fastq files after job done.
'''

data="/ltr/Fraietta/data"
out="/data2/count"
fastq="/home/ubuntu/data/viralLTR/fastq"
max_jobs=2
for dir in "16CT023-14_TDN" "19CT023-04_TDN"  \
"19CT023-12_TDN" "16CT023-12_D7" "19CT023-02_D7" "19CT023-09_D10"  \
"19CT023-18_D7" "16CT023-12_TDN" "19CT023-02_TDN" "19CT023-09_TDN" "19CT023-18_TDN"; do
    
    while [[ $(jobs -p | wc -l) -ge $max_jobs ]]; do
        sleep 10
    done
    
    (
		logFile=${out}/${dir}_log.txt
		touch $logFile
		echo ${data}/${dir} >> $logFile
		mkdir ${fastq}/${dir}
		cat ${data}/${dir}/R1* > ${fastq}/${dir}/${dir}"_S1_L001_R1_001.fastq.gz"
		cat ${data}/${dir}/R2* > ${fastq}/${dir}/${dir}"_S1_L001_R2_001.fastq.gz"
		echo "fastq files generated" >> $logFile
		cd ${out}
		cellranger count --id=$dir \
		--transcriptome=/home/ubuntu/tool/cellRange/referenceGenome/refdata-gex-GRCh38-2020-A \
		--fastqs=${fastq}/${dir} --localcores=16 --localmem=124 --no-bam --nosecondary 
		echo ${dir}" cellranger completed" >> $logFile
		rm -r ${fastq}/${dir}
		echo "removed fastq folder"  >> $logFile
	) &
done

wait
echo "counting completed"