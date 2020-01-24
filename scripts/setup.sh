#!/usr/bin/env bash
DIR="build"
[ ! -d "$DIR" ] && mkdir -p "$DIR"
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=./${CMTCONFIG}/${LCG_VERSION}/ && make -j12 && make install
