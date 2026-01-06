[WORKSPACE INSTALLATION INSTRUCTIONS]
Note from Sage: The PyOorb library is only available for Mac and Linux. To get around this we'll have to use a WSL environment. This readme.txt will help you set that up and get everything working on your local computer.

1) In the Windows Command Prompt, enter the following command to install WSL: wsl --install
2) It should prompt you for a UNIX username + pass. Enter whatever you'd like (I did my NAU username for clarity)
3) WSL should now be installed. You should have an application called Ubuntu
4) Boot up Ubuntu. It will be a terminal. Input all remaining commands into this terminal.
5) Run: wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
6) Run: bash Miniconda3-latest-Linux-x86_64.sh
7) Restart Ubuntu
8) Conda should now be installed.
9) Run "explorer.exe ." in the terminal to pull up your Linux file system location
	^ Note: You should save this location because this will be your working environment. I added a shortcut to it to my desktop.
10) Put the environments.yml file in this folder (or make an LSST folder and put it inside there)
11) Create a new Conda environment (this will be your LSST_Project environment)
12) Activate that Conda environment

13) RUN THE FOLLOWING COMMANDS (from your project folder):
conda config --add channels conda-forge
conda install numpy
conda install matplotlib
conda install spiceypy
conda install pyarrow
conda install pandas
conda install astropy
conda install sbpy
conda install pymongo
conda install openorb

14) Restart Ubuntu

15) You should now be able to run main.py 	(but make sure you reactivate your conda environment and CD into the workspace first)

16) From there, every time you boot up Ubuntu you need to:
	a) Activate your Conda environment
	b) cd into your project folder
	c) Run whatever .py files you'd like

