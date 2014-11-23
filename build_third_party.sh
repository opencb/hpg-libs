#!/bin/bash

git submodule init
git submodule update
cd third_party/htslib/
make
cd ../samtools
make
cd ../..

