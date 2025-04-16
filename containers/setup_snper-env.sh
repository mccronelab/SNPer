#! bin/bash

# create venv
python3 -m venv snper-env \
# activate it
source snper-env/bin/activate

pip installl --upgrade pip
pip install setuptools
# install Bio for convert_tsv_coords.py
pip install Bio

# install liftoff
git clone https://github.com/agshumate/Liftoff liftoff 
cd liftoff
python setup.py install

# install cutadapt
pip install cutadapt