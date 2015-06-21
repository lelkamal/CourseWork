#!/bin/bash
#SBATCH --output=p1.out
#SBATCH --job-name=p1
#SBATCH --error=p1.err
#SBATCH --nodes=1

MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

mpiexec.hydra --bootstrap slurm -f $MACHINEFILE  -env I_MPI_DEVICE rdma:ofa-v2-mthca0-1u -ppn 1 -n 1 ./a.out

rm $MACHINEFILE

