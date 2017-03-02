#!/bin/bash

#PBS -N me
#PBS -o job_reports/me-output
#PBS -e job_reports/me-error
#PBS -l walltime=12:00:00
#PBS -t 1-16
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
	echo "${1}"
	PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
splits=16

Rscript \
	mz_variance.r \
	${i} \
	${splits} \
	../results/mz${i}.rdata
