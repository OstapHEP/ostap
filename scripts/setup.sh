#!/usr/bin/env bash
DIR="build"
[ ! -d "$DIR" ] && mkdir -p "$DIR"
cd $DIR
cmake .. -DCMAKE_INSTALL_PREFIX=./INSTALL/ -GNinja && ninja && ninja install
