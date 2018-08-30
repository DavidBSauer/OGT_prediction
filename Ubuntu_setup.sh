#!/bin/bash

#install fresh ubuntu

#install tRNAscan-SE 
echo "installing tRNAscan"
sudo apt-get update
sudo apt-get install gcc make perl automake
wget http://eddylab.org/infernal/infernal-1.1.2.tar.gz
tar -zxvf infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure
make
make check
sudo make install
cd ..

wget http://trna.ucsc.edu/software/trnascan-se-2.0.0.tar.gz
tar -zxvf trnascan-se-2.0.0.tar.gz
cd tRNAscan-SE-2.0
./configure
make
sudo make install
sudo chmod 755 /usr/local/lib/tRNAscan-SE/tRNAscanSE/
cd ..

#install barrnap-0.9
echo "installing barrnap"
cd $HOME
wget https://github.com/tseemann/barrnap/archive/0.9.tar.gz
tar -zxvf 0.9.tar.gz
echo "PATH=$PATH:$HOME/barrnap-0.9/bin" >> .bashrc

#install python
echo "installing python and necessary python packages"
sudo apt-get -y install python-pip python-dev build-essential
sudo apt-get -y install bedtools

#install python packages
sudo pip install numpy scipy matplotlib biopython bcbio-gff tqdm sklearn matplotlib-venn

#install GeneMarkS
#go to exon.gatech.edu/GeneMark/license_download.cgi
#download Genemark-S and the 64_bit key
echo "Don't forget:"
echo "1. download and install GeneMarkS"
echo "2. edit the external_tools.txt file as needed"
