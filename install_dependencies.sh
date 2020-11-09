#!/bin/bash

git clone https://github.com/egonozer/AGEnt --branch v0.3.1
git clone https://git.code.sf.net/p/clustage/code --branch v0.9 ClustAGE

wget http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar -zxvf minikraken_20171019_8GB.tgz
mv -r minikraken_20171019_8GB ref/

#wget -O ref/UniVec https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
#makeblastdb -in ref/UniVec -dbtype nucl
