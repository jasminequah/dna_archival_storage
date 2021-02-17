conda create -n dna_project python pip pytorch=1.5.0 torchvision cudatoolkit 'numpy<=1.18.5' -c pytorch -c conda-forge
conda activate dna_project
pip3 install scrappie fast5_research
cd bonito
pip3 install -r requirements.txt
python3 setup.py develop
bonito download --ctc
