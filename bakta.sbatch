#!/bin/bash

#SBATCH -J bakta
#SBATCH -p cpu_p
#SBATCH -q cpu_short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --mail-user=roberto.siani
#SBATCH --mail-type=ALL
#SBATCH --nice=10000

source $HOME/.bashrc

conda activate bakta

echo "Started bakta annotation"

bakta blic.fna --db ~/db --output blic_bakta  --prefix blic --genus Bacillus --species licheniformis --strain AK1 --gram "?" --locus bl --locus-tag bl  -t 8 
bakta bvel.fna --db ~/db --output bvel_bakta  --prefix bvel --genus Bacillus --species velezensis --strain FZB42 --gram "?"  --locus bv --locus-tag bv -t 8

