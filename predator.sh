#!/bin/bash
#SBATCH --output=Predator.out
#SBATCH --job-name=Predator
#SBATCH --error=Predator.err
#SBATCH --nodes=1

MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

mpiexec.hydra --bootstrap slurm -f $MACHINEFILE  -env I_MPI_DEVICE rdma:ofa-v2-mthca0-1u -ppn 4 -n 4 ./a.out

rm $MACHINEFILE

