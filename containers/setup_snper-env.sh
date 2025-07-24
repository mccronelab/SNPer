#! bin/bash

# create venv
python3 -m venv snper-env \
# activate it
source snper-env/bin/activate

pip install --upgrade pip
pip install setuptools
# install Bio for convert_tsv_coords.py
pip install Bio

# install liftoff
git clone https://github.com/agshumate/Liftoff liftoff 
cd liftoff
pip install .
cd

# install cutadapt
pip install cutadapt