#!/bin/bash

APP_NAME=lammps
DOCDIR="$1"
VERSION="$2"
BRANCH=$(git branch --show-current)
DESTDIR=${PWD}

echo "Delete tar files, if they exist"
rm -rvf ${PWD}/lammps-src-*.tar.gz
pushd ${DOCDIR}/..
git archive --output=${DESTDIR}/${BRANCH}.tar --prefix=lammps-${VERSION}/ HEAD

cd ${DOCDIR}
make clean-all
make html pdf
tar -rf ${DESTDIR}/${BRANCH}.tar --transform "s,^,lammps-${VERSION}/doc/," html Manual.pdf
popd
mv ${BRANCH}.tar lammps-src-${VERSION}.tar
gzip -9v lammps-src-${VERSION}.tar
exit 0
