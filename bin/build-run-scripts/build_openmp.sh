#!/usr/bin/env bash
module load openmpi/4.1.2

VERSION="imdv3"

PROJECT_ROOT=$(pwd)/../..
BUILD_DIR="$PROJECT_ROOT/bin/build-openmp"
RUN_DIR="$PROJECT_ROOT/bin/run-openmp"

IMD_ASYNC="ON"            # ON or OFF

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR" "$RUN_DIR"

EXEC_NAME="lmp_openmp_${VERSION}_${IMD_ASYNC}"

cd "$PROJECT_ROOT"
git checkout "${VERSION}-benchmarking"

# LAMMPS with MPI + OpenMP only (no CUDA/Kokkos/GPU)
cmake -S "$PROJECT_ROOT/cmake" -B "$BUILD_DIR" \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS_RELEASE="-march=native" \
  -D CMAKE_C_FLAGS_RELEASE="-march=native" \
  -D BUILD_MPI=OFF \
  -D BUILD_OMP=ON \
  -D PKG_GPU=OFF \
  -D PKG_CHARMM=ON \
  -D PKG_CMAP=ON \
  -D PKG_COLVARS=ON \
  -D PKG_SHAKE=ON \
  -D PKG_KSPACE=ON \
  -D PKG_RIGID=ON \
  -D PKG_MOLECULE=ON \
  -D PKG_MISC=ON \
  -D LAMMPS_ASYNC_IMD=${IMD_ASYNC} \
  -D PKG_OPENMP=ON \
  -D PKG_EXTRA-FIX=ON \
  -D PKG_EXTRA-DUMP=ON \
  -D PKG_MOLFILE=ON \
  -D LAMMPS_MACHINE=openmp

cmake --build "$BUILD_DIR" -- -j16

# copy executable to run dir
cp "$BUILD_DIR/lmp_openmp" "$RUN_DIR/$EXEC_NAME"