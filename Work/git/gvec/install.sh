cd /tokp/work/jpuch/git/gvec

module purge
module load git/2.43
module load git-lfs/3.3
module load cmake/3.30
module load intel/2023.1.0.x
module load impi/2021.9
module load impi-interactive/1.0
module load mkl/2023.1
module load hdf5-serial/1.14.1
module load netcdf-serial/4.9.2
module load python-waterboa/2024.06
module list
export FC=`which mpiifort`
export CC=`which mpiicc`
export CXX=`which mpiicc`

git clone git@gitlab.mpcdf.mpg.de:gvec-group/gvec.git gvec
cd gvec
python -m venv .venv
source .venv/bin/activate
pip install .[dev,examples] -v
