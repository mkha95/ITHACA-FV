#!/bin/bash
source /usr/lib/openfoam/openfoam2106/etc/bashrc && \
source etc/bashrc && \
./Allwmake -j 4 && \
apt-get update && \
apt-get install -y cmake git && \
git clone https://github.com/google/googletest.git && \
cd googletest && \
mkdir build           && \
cd build && \
cmake ..              && \
make && \
sudo make install &&
export GTEST_DIR='/usr/local/lib' &&  \
cd ../.. && \
wmake  && \
./test_all.exe ;

