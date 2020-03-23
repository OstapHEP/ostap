#!/usr/bin/env bash
DIR="LCG_${LCG_VERSION}"
[ ! -d "$DIR" ] && mkdir -p "$DIR"
cd $DIR
cmake .. -DCMAKE_INSTALL_PREFIX=./INSTALL/ && make -j12 && make install
