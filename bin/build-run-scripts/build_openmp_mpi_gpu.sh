#!/usr/bin/env bash
module load openmpi/4.1.2
module load cuda/12.3/ubuntu22.04/64

PROJECT_ROOT=$(pwd)/../..
BUILD_DIR="$PROJECT_ROOT/bin/build-openmp-mpi-gpu"
RUN_DIR="$PROJECT_ROOT/bin/run-openmp-mpi-gpu"
EXEC_NAME="lmp_openmp_mpi_gpu"
CUDA_ROOT="/nfs/packages/opt/Linux_x86_64/cudatoolkit/12.3/ubuntu22.04"

# GPU Architecture selection
#GPU_ARCH="sm_80"         # NVIDIA A100 80GiB
#GPU_ARCH="sm_86"         # NVIDIA A30 24GiB  
GPU_ARCH="sm_89"          # NVIDIA RTX 4070 Ti (Ada Lovelace)
#GPU_ARCH="sm_90"         # NVIDIA GH200 (Grace Hopper)

# Kokkos Architecture
#KOKKOS_ARCH="AMPERE80"   # For A100 (sm_80)
#KOKKOS_ARCH="AMPERE86"   # For A30 (sm_86)
KOKKOS_ARCH="ADA89"       # For RTX 4070 Ti (sm_89)
#KOKKOS_ARCH="HOPPER90"   # For GH200 (sm_90)

GPU_PREC="single"          # double or single or mixed precision

rm -rf "$BUILD_DIR"

mkdir -p "$BUILD_DIR" "$RUN_DIR"

# LAMMPS with OpenMP, MPI, GPU and Kokkos

# removed -ffast-math \
cmake -S "$PROJECT_ROOT/cmake" -B "$BUILD_DIR" \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS_RELEASE="-march=native" \
  -D CMAKE_C_FLAGS_RELEASE="-march=native" \
  -D BUILD_MPI=ON \
  -D BUILD_OMP=ON \
  -D PKG_GPU=ON \
  -D PKG_KOKKOS=ON \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ARCH_${KOKKOS_ARCH}=ON \
  -D CMAKE_CXX_COMPILER="$PROJECT_ROOT/lib/kokkos/bin/nvcc_wrapper" \
  -D FFT_KOKKOS=CUFFT \
  -D Kokkos_ENABLE_DEPRECATION_WARNINGS=OFF \
  -D PKG_CHARMM=ON \
  -D PKG_CMAP=ON \
  -D PKG_SHAKE=ON \
  -D PKG_KSPACE=ON \
  -D PKG_RIGID=ON \
  -D PKG_MOLECULE=ON \
  -D PKG_OPENMP=ON \
  -D LAMMPS_MACHINE=openmp_mpi_gpu \
  -D GPU_API=cuda \
  -D GPU_ARCH=${GPU_ARCH} \
  -D FFT_GPU=CUFFT \
  -D GPU_PREC=$GPU_PREC \
  -D CUDA_TOOLKIT_ROOT_DIR="$CUDA_ROOT" \
  -D CMAKE_CUDA_COMPILER="$CUDA_ROOT/bin/nvcc"

cmake --build "$BUILD_DIR" -- -j12

# copy executable to run dir
cp "$BUILD_DIR/lmp_openmp_mpi_gpu" "$RUN_DIR/$EXEC_NAME"