#!/bin/bash

#PBS -N sibp
#PBS -o job_reports/sibp-output
#PBS -e job_reports/sibp-error
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
	sib_mr_popstrat.r \
	${i} \
	${splits} \
	../results/sib_popstrat_${i}.rdata
