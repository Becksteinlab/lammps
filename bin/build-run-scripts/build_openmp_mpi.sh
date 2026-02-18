#!/usr/bin/env bash
module load openmpi/4.1.2
# module load cmake-3.26.5-gcc-11.2.0

PROJECT_ROOT=$(pwd)/../..
BUILD_DIR="$PROJECT_ROOT/bin/build-openmp-mpi"
RUN_DIR="$PROJECT_ROOT/bin/run-openmp-mpi"
EXEC_NAME="lmp_openmp_mpi"

mkdir -p "$BUILD_DIR" "$RUN_DIR"

cmake -S "$PROJECT_ROOT/cmake" -B "$BUILD_DIR" \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS_RELEASE="-march=native" \
  -D CMAKE_C_FLAGS_RELEASE="-march=native" \
  -D BUILD_MPI=ON \
  -D BUILD_OMP=ON \
  -D PKG_GPU=OFF \
  -D PKG_CHARMM=ON \
  -D PKG_CMAP=ON \
  -D PKG_SHAKE=ON \
  -D PKG_KSPACE=ON \
  -D PKG_RIGID=ON \
  -D PKG_MOLECULE=ON \
  -D PKG_OPENMP=ON \
  -D LAMMPS_MACHINE=openmp_mpi

cmake --build "$BUILD_DIR" -- -j6

# copy executable to run dir
cp "$BUILD_DIR/lmp_openmp_mpi" "$RUN_DIR/$EXEC_NAME"

