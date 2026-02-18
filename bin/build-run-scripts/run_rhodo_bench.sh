#!/usr/bin/env bash

module load openmpi/4.1.2
module load cuda/12.3/ubuntu22.04/64

PROJECT_ROOT=$(pwd)/../..
BENCH_DIR="$PROJECT_ROOT/bench"
INPUT_FILE="$BENCH_DIR/in.rhodo"

max_cores=8
nompthreads=4 # Number of OpenMP threads from 1 to 8
nmpitasks=1 # Number of MPI tasks 1 to 8
ngpu=1 # only 1
nompthreads=$(($nompthreads > $max_cores/$nmpitasks ? $max_cores/$nmpitasks : $nompthreads)) # minimum of 8/ntasks or nompthreads

# Only bind when using single MPI task
if [ $nmpitasks -eq 1 ]; then
    export OMP_PROC_BIND=spread # close
else
    export OMP_PROC_BIND=false
fi

export OMP_PLACES=threads
export OMP_NUM_THREADS=$nompthreads

# OpenMP
RUN_DIR="$PROJECT_ROOT/bin/run-openmp"
EXEC_NAME="lmp_openmp"

mkdir -p "$RUN_DIR"

cd "$RUN_DIR" || exit 1

"./$EXEC_NAME" -sf omp -pk omp $nompthreads -in "$INPUT_FILE" -log "log.rhodo.omp${nompthreads}" > "log.rhodo.omp${nompthreads}.full" 2>&1

# OpenMP+MPI run
RUN_DIR="$PROJECT_ROOT/bin/run-openmp-mpi"
EXEC_NAME="lmp_openmp_mpi"

mkdir -p "$RUN_DIR"

cd "$RUN_DIR" || exit 1

mpirun -np $nmpitasks --bind-to core --map-by slot:PE=$nompthreads -x OMP_NUM_THREADS=$nompthreads \
            "./$EXEC_NAME" -sf omp -pk omp $nompthreads -in "$INPUT_FILE" -log "log.rhodo.mpi${nmpitasks}xomp${nompthreads}" > "log.rhodo.omp_mpi${nmpitasks}xomp${nompthreads}.full" 2>&1

# OpenMP+MPI+GPU run

RUN_DIR="$PROJECT_ROOT/bin/run-openmp-mpi-gpu"
EXEC_NAME="lmp_openmp_mpi_gpu_single"

mkdir -p "$RUN_DIR"

cd "$RUN_DIR" || exit 1

# Kokkoks run

mpirun -np $nmpitasks ./"$EXEC_NAME" -k on g 1 t $nompthreads -sf kk -pk kokkos newton on neigh half -in "$INPUT_FILE" -log log.rhodo.kokkos.mpi${nmpitasks}xomp${nompthreads}

# GPU run

mpirun -np $nmpitasks ./"$EXEC_NAME" -sf gpu -pk gpu 1 neigh yes newton off -in "$INPUT_FILE" -log log.rhodo.gpu.mpi${nmpitasks}xomp${nompthreads}

# GPU + OMP hybrid run

mpirun -np $nmpitasks ./"$EXEC_NAME" -sf hybrid gpu omp -pk gpu 1 neigh yes newton off -pk omp $nompthreads neigh no -in "$INPUT_FILE" -log log.rhodo.hybrid-gpu-omp.mpi${nmpitasks}xomp${nompthreads}
# mpirun -np $nmpitasks ./"$EXEC_NAME" -sf gpu -pk gpu 1 neigh yes newton off -pk omp $nompthreads neigh no -in "$INPUT_FILE" -log log.rhodo.hybrid-gpu-omp.mpi${nmpitasks}xomp${nompthreads}

# Kokkos + GPU hybrid run

mpirun -np $nmpitasks ./"$EXEC_NAME" -k on g 1 t $nompthreads -sf hybrid kk gpu -pk kokkos newton on neigh half -pk gpu 1 -in "$INPUT_FILE" -log log.rhodo.hybrid-kokkos-gpu.mpi${nmpitasks}xomp${nompthreads}