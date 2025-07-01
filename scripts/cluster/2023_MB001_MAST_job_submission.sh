#!/bin/bash
#SBATCH -N 1                        # number of nodes
#SBATCH -n 2                        # number of cores
#SBATCH --mem 256G                   # memory pool for all cores
#SBATCH -t 0-12:00:00               # runtime limit (D-HH:MM:SS)
#SBATCH -o /g/bernabeu-hd/batzilla/HPC_log/slurm.%N.%j.out          # STDOUT
#SBATCH -e /g/bernabeu-hd/batzilla/HPC_log/slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=alina.batzilla@embl.de # send-to address

###### Path to Local Libraies
R_LocalLibraryPath=/g/huber/users/batzilla/R/x86_64-pc-linux-gnu-library/4.1

RscriptPath=/g/bernabeu-hd/batzilla/scRNAseq_final

###### Working direcotory
wd=/g/bernabeu-hd/batzilla/scRNAseq_final

###### Name of Log file
logfile=2023_MB001_MAST_log.txt

############################################################################################################


if [ ! -e ${wd} ]; then
	mkdir ${wd}
fi
cd ${wd}

module load R/4.1.2-foss-2021b
R --version | tee ${wd}/${logfile}
R --vanilla --slave -e ".libPaths(c(\"${R_LocalLibraryPath}\",.libPaths()))" | tee -a ${wd}/${logfile}
R --slave -e "setwd(\"${wd}\")" | tee -a ${wd}/${logfile}

Rscript --slave ${RscriptPath}/2023_MB001_MAST.R | tee -a ${wd}/${logfile}