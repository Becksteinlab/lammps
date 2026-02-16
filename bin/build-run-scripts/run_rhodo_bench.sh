#!/usr/bin/env bash

module load openmpi/4.1.2
module load cuda/12.3/ubuntu22.04/64

PROJECT_ROOT=$(pwd)/../..
RUN_DIR="$PROJECT_ROOT/bin/run-openmp-mpi-gpu"
EXEC_NAME="lmp_openmp_mpi_gpu_single"
BENCH_DIR="$PROJECT_ROOT/bench"
INPUT_FILE="$BENCH_DIR/in.rhodo"

mkdir -p "$RUN_DIR"

cd "$RUN_DIR" || exit 1

nompthreads=4
nmpitasks=1
ngpu=1

# Only bind when using single MPI task
if [ $nmpitasks -eq 1 ]; then
    export OMP_PROC_BIND=spread # close
else
    export OMP_PROC_BIND=false
fi

export OMP_PLACES=threads
export OMP_NUM_THREADS=$nompthreads

# Kokkoks run

mpirun -np $nmpitasks ./"$EXEC_NAME" -k on g 1 t $nompthreads -sf kk -pk kokkos newton on neigh half -in "$INPUT_FILE" -log log.rhodo.kokkos.mpi${nmpitasks}xomp${nompthreads}

# GPU run

mpirun -np $nmpitasks ./"$EXEC_NAME" -sf gpu -pk gpu 1 neigh yes newton off -in "$INPUT_FILE" -log log.rhodo.gpu.mpi${nmpitasks}xomp${nompthreads}

# GPU + OMP hybrid run

mpirun -np $nmpitasks ./"$EXEC_NAME" -sf hybrid gpu omp -pk gpu 1 neigh yes newton off -pk omp $nompthreads neigh no -in "$INPUT_FILE" -log log.rhodo.hybrid-gpu-omp.mpi${nmpitasks}xomp${nompthreads}
# mpirun -np $nmpitasks ./"$EXEC_NAME" -sf gpu -pk gpu 1 neigh yes newton off -pk omp $nompthreads neigh no -in "$INPUT_FILE" -log log.rhodo.hybrid-gpu-omp.mpi${nmpitasks}xomp${nompthreads}

# Kokkos + GPU hybrid run

mpirun -np $nmpitasks ./"$EXEC_NAME" -k on g 1 t $nompthreads -sf hybrid kk gpu -pk kokkos newton on neigh half -pk gpu 1 -in "$INPUT_FILE" -log log.rhodo.hybrid-kokkos-gpu.mpi${nmpitasks}xomp${nompthreads}