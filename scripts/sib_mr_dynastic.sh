#!/bin/bash

#PBS -N sibd
#PBS -o job_reports/sibd-output
#PBS -e job_reports/sibd-error
#PBS -l walltime=12:00:00
#PBS -t 1-100
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
	echo "${1}"
	PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
splits=100

cd ${HOME}/repo/mrtwin_power/scripts

Rscript \
	sib_mr_dynastic.r \
	${i} \
	${splits} \
	../results/sib_dynastic_${i}.rdata
