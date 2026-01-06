A) In the Windows Command Prompt, enter the following command to install WSL: wsl --install
B) It should prompt you for a UNIX username + pass. Enter whatever you'd like (I did my NAU username for clarity)
C) WSL should now be installed. You should have an application called Ubuntu
D) Boot up Ubuntu. It will be a terminal. Input all remaining commands into this terminal.
E) Run: wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
F) Run: bash Miniconda3-latest-Linux-x86_64.sh
G) Restart Ubuntu
H) Conda should now be installed.
I) Run "explorer.exe ." in the terminal to pull up your Linux file system location
J) Put the environments.yml file in this folder (or make an LSST folder and put it inside there)
J) Create a new Conda environment (this will be your LSST_Project environment)
K) Activate that Conda environment
L) 

1) Install Anaconda Distribution via https://www.anaconda.com/download
2) Run Anaconda Prompt in the Windows search bar
3) CD to your project folder (where environments.yml is)
4) Create a new Conda environment (this will be your LSST_Project environment)
5) Activate that Conda environment

6) RUN THE FOLLOWING COMMANDS (from your project folder):
pip install pyorb # if this doesn't work try "pip install pyorb"

conda config --add channels conda-forge
conda install numpy
conda install matplotlib
conda install spiceypy
conda install pyarrow
conda install pandas
conda install astropy
conda install sbpy
conda install pymongo

7) If everything works you should be able to run main.py

8) From there, every time you boot up the program you should run something that looks like this from Anaconda Prompt:
conda activate lsst_env                      # conda activate lsst2
cd C:\Users\grays\Desktop\LSSTProject
python main.py