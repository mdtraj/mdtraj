apt-get install python-software-properties git -y
apt-add-repository ppa:fkrull/deadsnakes -y
apt-get update
apt-get install python-pip libfreetype6-dev libpng12-dev libhdf5-serial-dev \
    g++ libatlas-base-dev gfortran netcdf-bin libnetcdf-dev libffi-dev -y
