# setup for detsim
export CAL_DB=/mu2e/data/users/gandr/twist/caldb_ascii

##setup mu2e
source /grid/fermiapp/products/mu2e/setupmu2e-art.sh

setup -B gsl v1_15a
setup -B log4cpp v1_0
setup -B boost v1_53_0 -q +e2:+prof
setup -B root v5_34_05 -q +e2:+prof
