# setup for  mu2egpvm machines
export CAL_DB=/mu2e/data/users/gandr/twist/caldb_ascii

setup mu2e
###source /grid/fermiapp/products/mu2e/setupmu2e-art.sh

setup -B boost v1_57_0 -q +e7:+prof
setup -B root v5_34_25 -q +e7:+prof

#setup -B gsl v1_15a
#setup -B log4cpp v1_0
PATH=/mu2e/app/users/gandr/local/install-sl6-e7/bin:$PATH
LD_LIBRARY_PATH=/mu2e/app/users/gandr/local/install-sl6-e7/lib:$LD_LIBRARY_PATH
