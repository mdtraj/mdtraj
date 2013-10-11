apt-get install python-software-properties git -y
apt-add-repository ppa:fkrull/deadsnakes -y
apt-get update

apt-get install libfreetype6-dev libpng12-dev libhdf5-serial-dev \
    g++ libatlas-base-dev gfortran netcdf-bin libnetcdf-dev libffi-dev -y
echo "apt-get install ${python}-dev"
echo "apt-get install ${PYTHON}-dev"
apt-get install python2.7-dev python2.7-dev python3.3-dev
