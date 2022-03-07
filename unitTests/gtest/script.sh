#!/bin/bash

source /usr/lib/openfoam/openfoam2106/etc/bashrc && \
    source etc/bashrc && \
    ./Allwmake -au -j 4;

