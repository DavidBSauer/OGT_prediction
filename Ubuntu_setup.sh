#!/bin/bash

#install fresh ubuntu

#install tRNAscan-SE 
echo "installing tRNAscan"
#sudo apt-add-repository -y "deb http://http.us.debian.org/debian sid main non-free"
sudo apt-get -y update
sudo apt-get -y install trnascan-se

#install barrnap-0.8
echo "installing barrnap"
cd $HOME
wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz
tar zxvf 0.8.tar.gz
echo "PATH=$PATH:$HOME/barrnap-0.8/bin" >> .bashrc

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
