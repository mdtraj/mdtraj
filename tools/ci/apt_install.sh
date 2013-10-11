sudo apt-get install python-software-properties git -y
sudo apt-add-repository ppa:fkrull/deadsnakes -y
sudo apt-get update

sudo apt-get install libfreetype6-dev libpng12-dev libhdf5-serial-dev \
    g++ libatlas-base-dev gfortran netcdf-bin libnetcdf-dev libffi-dev -y
sudo apt-get install ${python}-dev

