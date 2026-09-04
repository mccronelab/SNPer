#! bin/bash

# create venv
python3 -m venv snper-env \
# activate it
source snper-env/bin/activate

pip install --upgrade pip
pip install setuptools
# install Bio for the bin/ scripts
pip install Bio

# install cutadapt
pip install cutadapt