#!/bin/bash
INSTALL_DIR=/path/to/installation_location/ # CHANGE THIS LINE!!!

cd $INSTALL_DIR
TMP_DIR=$INSTALL_DIR/TMP/
mkdir $TMP_DIR
git clone https://github.com/Kalan-Lab/lsaBGC
conda env create -f $INSTALL_DIR/lsaBGC/lsaBGC_env.yml -p $INSTALL_DIR/lsaBGC_conda_env/
conda activate $INSTALL_DIR/lsaBGC_conda_env
cd $INSTALL_DIR/lsaBGC/
python setup.py install
pip install .
setup_annotation_dbs.py
echo $INSTALL_DIR
conda deactivate
echo $'To activate conda environment and use lsaBGC in the future, simply type:\nconda activate '"$INSTALL_DIR"'/lsaBGC_conda_env/'