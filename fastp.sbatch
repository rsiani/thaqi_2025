#!/bin/bash
#
#SBATCH --job-name=trimming
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=roberto.siani
#SBATCH --partition=cpu_p
#SBATCH --qos cpu_normal
#SBATCH --time=48:00:00
#SBATCH --nice=10000

source $HOME/.bashrc

conda activate fastp

alias java='java -d64'

echo "FASTP: trimming and quality filtering"

mkdir fastp_$(date +"%d%m")

parallel -j 10 -i {} --tempdir /localscratch/ --compress fastp --in1 {1} --in2 {2} --out1 fastp_$(date +"%d%m")/{1} --out2 fastp_$(date +"%d%m")/{2} -j fastp_$(date +"%d%m")/{1.}.json -h fastp_$(date +"%d%m")/{1.}.html --correction -p -t 1 --dedup -g -x -y -r ::: *R1_001.fastq.gz :::+ *R2_001.fastq.gz

echo "All done!!!"

