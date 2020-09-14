#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

//[[Rcpp::depends(BH)]]

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <string>
#include "brent.h"
#include <random> // for distributions
#ifdef _OPENMP
#   include <omp.h>
#endif
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdexcept> // for handling exceptions
#include "boost/math/special_functions/gamma.hpp" // for the incomplete gamma function
#include "model.h"
#include "utils.h"

using namespace std;
using namespace brent; // for zero()

using boost::math::tgamma_lower;


// set up random number generators
extern mt19937 generator;

const double check_tolerance(1e-15);

double Teqn(double Ts, double params[]) {
	double A = params[0];
	double b = params[1];
	double C = params[2];
	double n = params[3];
	double s0 = params[4];
	double k = params[5];
	double mu = params[6];

	double As = A*k*Ts;
	double Cs = C*k*Ts;

	double intfac = 1/mu * pow(Cs, n)* Ts/(n+1) * pow(1-s0,n+1);
	double bnfr = (b+1)/(n+1);

	double result = pow(As,b) / pow(Cs, n*bnfr) * pow(mu/Ts * (n+1),(b-n)/(n+1)) * tgamma_lower(bnfr, intfac) - exp(-intfac);

	return result;
}

// compute the constants for Tc
vector<double> const_Tc(double Ts, double T0, double params[]){
	double A = params[0];
	double b = params[1];
	double C = params[2];
	double n = params[3];
	double s0 = params[4];
	double k = params[5];
	double mu = params[6];

	double As = A*k*Ts;
	double Cs = C*k*Ts;

	double intfac = 1/mu * pow(Cs, n)* Ts/(n+1) * pow(T0/Ts-s0,n+1);
	double bnfr = (b+1)/(n+1);

	vector<double> res(4);
	// the constants are C1, C2, C3 and H^star(T0)

	res[0] = 1/mu * pow(As*max(T0/Ts - s0, 0.0), b);
	res[1] = 1/mu * pow(Cs*max(T0/Ts - s0, 0.0), n);
	res[2] = exp(-res[1] * T0 + intfac) * pow(As,b) / pow(Cs, n*bnfr) * pow(mu/Ts * (n+1),(b-n)/(n+1)) * tgamma_lower(bnfr, intfac);
	res[3] = exp(-res[1] * T0);

	return res;
}

vector<double> CANL(vector<double> params, int N, double w0, double k, double eps) {
	// w0 is the constant load; if it is negative, it is a ramp load test
	double mu_A = params[0];
	double s_A = params[1];
	double mu_b = params[2];
	double s_b = params[3];
	double mu_C = params[4];
	double s_C = params[5];
  	double mu_n = params[6];
	double s_n = params[7];
	double mu_s0 = params[8];
	double s_s0 = params[9];

	double t = r8_epsilon ( );

	double mu_s = 1;
	double T0 = w0/k;

	const double uprbd = 2000000; //if the generate number is nan, the result is set to be this number (needed to be increased since 4 years is already larger than 20000)

	double tmp;

	vector<double> res(N);

	// Generate random variates
	std::normal_distribution<double> normA(mu_A, s_A);
	std::normal_distribution<double> normb(mu_b, s_b);
	normal_distribution<double> normC(mu_C, s_C);
	normal_distribution<double> normn(mu_n, s_n);
	normal_distribution<double> norms0(mu_s0, s_s0);

  #ifdef _OPENMP
	  omp_set_dynamic(0);
  #endif
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < N; i++) {
	    double tpar[7];
	    tpar[0] = exp(normA(generator));
	    tpar[1] = exp(normb(generator));
	    tpar[2] = exp(normC(generator));
	    tpar[3] = exp(normn(generator));
	    double tempx = exp(norms0(generator));
	    tpar[4] = tempx / (1.0+tempx);
		tpar[5] = k; // the loading rate
	    tpar[6] = mu_s; //mu_s can be set to an appropriate constant (not estimated).

		try{
			double z = zero ( 0.00001, 0.1, t, Teqn, tpar);

			if ((Teqn(0.00001,tpar) > 0) - (Teqn(0.1,tpar)>0) == 0) {
				res[i] = uprbd;  // no solution within usual timeframe (0.1 hours = 38840psi) so set to upbrd.  otherwise Tc can incorrectly be set to zero.
			} else {
				if(T0 < 0.0 || z < T0){
					res[i] = z;
				}
				else{
					if(T0/z < tpar[4]){
						res[i] = uprbd;
					}
					else{
						vector<double> C = const_Tc(z, T0, tpar);
						tmp = -1/C[1] * log((C[0]/C[1]*C[3] + C[2]) / (1 + C[0]/C[1]));
						if(std::isnan(tmp)) res[i] = uprbd;
						else res[i] = tmp;
					}
				}
			}
		}
		catch(const exception& e)
		{
//			cout << "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
//			cout << tpar[0] << " " << tpar[1] << " " << tpar[2] << " " << tpar[3] << " " << tpar[4] << " " << tpar[5] << " " << tpar[6] << endl;
		}
	}
	return res;
}

vector<double> summary_stat(vector<double> data, double censored, vector<double> probs){

	int n_s = probs.size();
	vector<double> s_obs(n_s, 0);
	vector<double> data_censored;
	double count_censored = 0;

	for(int i = 0; i < data.size(); ++i){
		if(data[i] > censored || std::isnan(data[i])) count_censored += 1;
		else{
			data_censored.push_back(log(data[i]));
		}
	}
	if(data_censored.size() != 0){
		s_obs = quantile(data_censored, probs);
	}
	return(s_obs);
}

