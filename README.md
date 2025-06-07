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
  -remove lines starting with "module"  
  -at the first line of both files, add:  
```
source /opt/intel/oneapi/setvars.sh  
make -f makefile.prep clean  
(make -f makefile clean)  
```
## 4. Edit makefile and makefile.prep:
  -mpiifort -> mpiifx
  -at the end of makefile.prep, add:
```
clean:
	rm -f $(prep_o) xprep
```
  -at the end of makefile, add:
```
clean:
	rm -f $(fdns_o) xfdns
```
