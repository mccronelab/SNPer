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

# install cutadapt
pip install cutadapt

RUN apt-get purge -y --auto-remove \
        build-essential \
        python3-dev libpython3.*-dev \
        libparasail-dev libsigsegv-dev \
        git && \
    rm -rf /var/lib/apt/lists/* /root/.cache/pip