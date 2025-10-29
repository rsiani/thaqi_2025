#!/bin/bash
#
#SBATCH --job-name=kallistoQuant
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=roberto.siani
#SBATCH --partition=cpu_p
#SBATCH --qos cpu_short
#SBATCH --time=05:00:00
#SBATCH --nice=10000

source $HOME/.bashrc
conda activate kallisto


echo "kallisto quantification"
mkdir kall31_$(date +"%d%m")

cat bvel.list | while read i; do kallisto quant -i ~/stefanie/bvel_k31.idx -o kall31_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat blic.list | while read i; do kallisto quant -i ~/stefanie/blic_k31.idx -o kall31_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done


echo "All done!!!"

