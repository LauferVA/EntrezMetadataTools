srun --ntasks=1 --cpus-per-task=1 --mem-per-cpu=65535 --partition=long --job-name=Entrez_Metadata --pty /bin/bash

module load Python/3.10.8-GCCcore-12.2.0

pip3 install Bio
pip3 install bs4 lxml html5lib
pip3 install certifi ssl
