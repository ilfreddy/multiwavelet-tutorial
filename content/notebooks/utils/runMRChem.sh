#!/bin/bash

# Set defaults
JOBNAME="jobname"
NTASKS="8"
NCPUS="8"
TIMELIMIT="10:00:00"
MEMORY="60"
SUBMIT=false

# Define usage functions
usage() { 
  echo "Usage: $0 [ -h ] [ -i JOBNAME ] [ -t NTASKS ] [ -c NCPUS ] [ -T TIMELIMIT ] [ -m MEMORY ] [ -x ]"
  echo
  echo "Optional arguments:"
  echo "-h                : Show this message and exit 0"
  echo "-i <JOBNAME>      : Input file w.o. extension (.inp assumed)  default: ${JOBNAME}"
  echo "-t <NTASKS>       : Number of MPI processes                   default: ${NTASKS}"
  echo "-c <NCPUS>        : Number of OMP threads                     default: ${NCPUS}"
  echo "-T <TIMELIMIT>    : Wall time limit in SLURM format           default: ${TIMELIMIT}"
  echo "-m <MEMORY>       : Total memory in GB                        default: ${MEMORY}"
  echo "-x                : Submit job to queue                       default: ${SUBMIT}"
}
bad_exit() { 
  usage
  exit 1 
}  
    
#    Parse arguments
while getopts "hi:t:c:T:m:x" options; do
  case "${options}" in
    h)
      usage
      exit 0;;
    i)
      JOBNAME=${OPTARG};;
    t)
      NTASKS=${OPTARG};;
    c)
      NCPUS=${OPTARG};;
    T)
      TIMELIMIT=${OPTARG};;
    m)
      MEMORY=${OPTARG};;
    x) 
      SUBMIT=true;;
    :)
      bad_exit;;
    *)
      bad_exit;;
  esac
done

# Check if file exists
if [ -f "${JOBNAME}.job" ]; then
  echo "File ${JOBNAME}.job already exists. Overwrite? ([y]/n)"
  read ANSWER
  if [ -z "${ANSWER}" ] || [ "${ANSWER}" == "y" ] || [ "${ANSWER}" == "yes" ]; then
    :
  else
    echo Aborted
    exit 1
  fi
fi

# Write SLURM job file
echo "#!/bin/bash" > ${JOBNAME}.job
echo "" >> ${JOBNAME}.job
echo "#SBATCH --account=nn4654k" >> ${JOBNAME}.job
echo "#SBATCH --mail-type=None" >> ${JOBNAME}.job
echo "#SBATCH --job-name=${JOBNAME}" >> ${JOBNAME}.job
echo "#SBATCH --output=${JOBNAME}.log" >> ${JOBNAME}.job
echo "#SBATCH --error=${JOBNAME}.err" >> ${JOBNAME}.job
echo "#SBATCH --time=${TIMELIMIT}" >> ${JOBNAME}.job
echo "#SBATCH --mem=${MEMORY}GB" >> ${JOBNAME}.job
echo "#SBATCH --ntasks=${NTASKS}" >> ${JOBNAME}.job
echo "#SBATCH --cpus-per-task=${NCPUS}" >> ${JOBNAME}.job
echo  >> ${JOBNAME}.job
echo "module purge" >> ${JOBNAME}.job
echo "module load MRChem/oieoifs" >> ${JOBNAME}.job
echo  >> ${JOBNAME}.job
echo "cp ${JOBNAME}.inp $SCRATCH" >> ${JOBNAME}.job
echo "cd $SCRATCH" >> ${JOBNAME}.job
echo "export OMP_NUM_THREADS=${NCPUS}" >> ${JOBNAME}.job
echo "mrchem --launcher 'srun -n ${NTASKS}' ${JOBNAME} > ${JOBNAME}.out" >> ${JOBNAME}.job
echo  >> ${JOBNAME}.job
echo "savefile ${JOBNAME}.out" >> ${JOBNAME}.job
echo "savefile ${JOBNAME}.json" >> ${JOBNAME}.job
echo "exit 0" >> ${JOBNAME}.job

echo "Job file written: ${JOBNAME}.job"

if [ "${SUBMIT}" = true ]; then
  sbatch ${JOBNAME}.job
fi

exit 0
