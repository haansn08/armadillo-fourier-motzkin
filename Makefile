arma_headers_location = ~/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include

ALL:
	g++ -Wall fourier-motzkin.cpp -o fourier-motzkin -g -I $(arma_headers_location) -DARMA_DONT_USE_WRAPPER -lblas -llapack
