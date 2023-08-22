#!/bin/bash

git clone https://github.com/egonozer/AGEnt --branch v0.3.1
git clone https://git.code.sf.net/p/clustage/code --branch v0.9 ClustAGE

wget http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar -zxvf minikraken_20171019_8GB.tgz
mv -r minikraken_20171019_8GB ref/

curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-adaptor.sif
mv fcs-adaptor.sif scripts/
