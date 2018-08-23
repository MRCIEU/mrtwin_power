#!/bin/bash


#SBATCH --job-name=dynastic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=201-300
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=mrcieu
#SBATCH --mem=4G

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

splits=300
set -e

mkdir -p ../results/sib_dynastic/
Rscript \
	sib_mr_dynastic.r \
	${i} \
	${splits} \
	../results/sib_dynastic/sib_dynastic_${i}.rdata
