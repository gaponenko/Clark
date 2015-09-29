# setup for detsim
export CAL_DB=/mu2e/data/users/gandr/twist/caldb_ascii

##setup mu2e
source /grid/fermiapp/products/mu2e/setupmu2e-art.sh

setup gsl v1_15a
setup log4cpp v1_0
setup boost v1_53_0 -q e2:prof
setup root v5_34_05 -q e2:prof
