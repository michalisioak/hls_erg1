#!/bin/bash
PROJECT_DIR=$(pwd)
BUILD_DIR="${PROJECT_DIR}/build"
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}
cmake ${PROJECT_DIR}
make
cp ${BUILD_DIR}/erg1 ${PROJECT_DIR}/
echo "Build complete! Executable has been copied to the project root."
