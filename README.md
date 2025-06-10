# Running fdnsrfv-mpi in Ubuntu:
## 1. Download intel oneAPI:
**download the key to system keyring**
```
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB -O key.pub  
gpg --dearmor < key.pub > key.gpg  
sudo cp key.gpg /usr/share/keyrings/oneapi-archive-keyring.gpg  
```
**add signed entry to apt sources and configure the APT client to use Intel repository:**
```
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
```
**update package list**
```
sudo apt update
```
**install hpc toolkit**
```
sudo apt install intel-hpckit
```
## 2. Install GNU-make
```
sudo apt install make
```
## 3. Edit build.sh and build-prep.sh:
  - remove lines starting with "module"  
  - at the first line of both files, add:  
```
source /opt/intel/oneapi/setvars.sh  
make -f makefile.prep clean  
(make -f makefile clean)  
```
## 4. Edit makefile and makefile.prep:
  - mpiifort -> mpiifx  
  - at the end of makefile.prep, add:
```
clean:
	rm -f $(prep_o) xprep
```
  - at the end of makefile, add:
```
clean:
	rm -f $(fdns_o) xfdns
```
## 5. Compile xprep and xfdns:  
  - enter the following:
```
chmod +x ./build-prep.sh
./build-prep.sh
chmod +x ./build.sh
./build.sh
```
## 6. Compile prefluid.ex and tecout.ex:  
```
source /opt/intel/oneapi/setvars.sh  
mpiifx -o ./prefluid.ex ./prefluid.f  
mpiifx -o ./tecout.ex ./tecout.f  
```
## 7. Run xprep to split input files:
```
source /opt/intel/oneapi/setvars.sh (enter this line when you open a new command prompt)
mpirun -n <number of processors> ./xprep
```
## 8. Run xfdns:
```
mpirun -n <number of processors> ./xfdns
```
# Running fdnsrfv-mpi in Windows:
## 1. Compile fdns:
  1. Download Visual Studio and Intel oneAPI
  2. Copy the following files to your fdns directory:
```
build.bat
build-prep.bat
makefile
makefile.prep
```
  3. Open Intel oneAPI command prompt before compiling
  4. Change to target directory
  5. Entering build-prep.bat (build.bat)
## 2. Run xprep and xfdns:
  1. Open Intel oneAPI command prompt
  2. entering:
```
mpiexec -n <number of processors> ./xprep.exe
(mpiexec -n <number of processors> ./xfdns.exe)
```
## Note:
  - shell -> batch
  - make -> nmake
# Run the main program:
  1. Check number of processor in **card.py**
  2. Edit the settings relative to the solver
  3. Execute **meshing.py, card.py, flowfield.py** to generate fort.11~13 (meshing.py needs to be run first)
  4. Execute **xprep(xprep.exe)** and split fort.11~13 in order
  5. Execute **xfdns(xfdns.exe)**
# Some Notes:
  1. IIQMAX in **fdns01** and **xprep.f** may need to be changed
  2. The size of IZS, IZT...... in **tecout.f** may need to be changed
  3. The direction of grids (i,j,k) needs to obey the **right hand rule**
  4. binary_fort13 in plot3dout.py, flowfield.py, tecout.py and test.py are of no use