#! bin/bash

# create venv
python3 -m venv liftoff-env \
# activate it
source liftoff-env/bin/activate

pip installl --upgrade pip
pip install setuptools

# install liftoff
git clone https://github.com/agshumate/Liftoff liftoff 
cd liftoff
python setup.py install