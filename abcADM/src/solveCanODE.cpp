#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

//[[Rcpp::depends(BH)]]

#include <stdlib.h>
#include <math.h>
#include <iostream> // for cin and Rcout
//#include <fstream>
//#include <sstream>
#include <vector> // for vector
#include <map>
#include <string>
#ifdef _OPENMP
#   include <omp.h>
#endif
#include <random> // for distributions
#include "boost/math/special_functions/gamma.hpp" // for the incomplete gamma function
#include "boost/numeric/odeint.hpp" // for the ode solver
#include "boost/math/tools/roots.hpp" // for the root finding
#include "utils.h"
#include <Rcpp.h>

using namespace std;
using namespace boost::numeric::odeint;
using boost::math::tgamma_lower;
using namespace boost::math::tools; // for bisect
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////

// setup parameters
const double mu = 1;
const double load_level = 4500;
const double k = 388440;
double T0 = load_level/k;
int n_per_theta = 1000;
const eps_tolerance<double> tol(6);
const double t_start = 0;
double t_end = 262980; // hours in 30 years
double constdt = 100;



// set up random number generator2
mt19937 generator2(0);

////////////////////////////////////////////////////////////////////////////////

double stepfun(double x, vector<double> x_vec, vector<double> y_vec){
	if(y_vec.size() - x_vec.size() != 1){
		Rcout << "The length of y_vec must be greater than the length of x_vec by 1." << endl;
		return(0);
	}
	int n = x_vec.size();
	if(x < x_vec[0]) return y_vec[0];
	else if(x >= x_vec[n-1]) return y_vec[n];
	else{
		for(int i = 1; i < n; ++i){
			if(x < x_vec[i] && x >= x_vec[i-1]) return y_vec[i];
		}
	}

	Rcout << "error: could not compute load profile value. \n";
	return 0;
}

class load_profile {

	double phi, load_s_shape, load_s_scale, load_p_shape, load_p_scale,  R0, alpha_d, alpha_l,
	load_d_mean, load_d_sd, gamma, mean_Ts, mean_Te, mean_Tp, D_d;
	int N = 250;


public:
	vector<double> T_s = vector<double>(N);
	vector<double> T_e = vector<double>(N);
	vector<double> load_s = vector<double>(N+1);
	vector<double> load_e = vector<double>(N+1);

	load_profile(map<string, double> load_profile_par):
	phi(load_profile_par["phi"]),
	load_s_shape(load_profile_par["load_s_shape"]),
	load_s_scale(load_profile_par["load_s_scale"]),
	load_p_shape(load_profile_par["load_p_shape"]),
	load_p_scale(load_profile_par["load_p_scale"]),
	R0(load_profile_par["R0"]),
	alpha_d(load_profile_par["alpha_d"]),
	alpha_l(load_profile_par["alpha_l"]),
	load_d_mean(load_profile_par["load_d_mean"]),
	load_d_sd(load_profile_par["load_d_sd"]),
	gamma(load_profile_par["gamma"]),
	mean_Ts(load_profile_par["mean_Ts"]),
	mean_Te(load_profile_par["mean_Te"]),
	mean_Tp(load_profile_par["mean_Tp"])
	{
		normal_distribution<double> norm_d(load_d_mean, load_d_sd);
		gamma_distribution<double> gamma_s(load_s_shape, load_s_scale);
		gamma_distribution<double> gamma_p(load_p_shape, load_p_scale);
		exponential_distribution<double> exp_Ts(1.0/mean_Ts);
		exponential_distribution<double> exp_Te(1.0/mean_Te);
		exponential_distribution<double> exp_Tp(1.0/mean_Tp);
		D_d = norm_d(generator2);
		T_s[0] = exp_Ts(generator2);
		T_e[0] = exp_Te(generator2);
		load_s[0] = gamma_s(generator2);
		load_e[0] = 0;
		for(int i = 1; i < N; ++i){
			T_s[i] = T_s[i-1] + exp_Ts(generator2); // times of owner changing
			load_s[i] = gamma_s(generator2);
			if(i % 2 != 0){
				T_e[i] = T_e[i-1] + exp_Tp(generator2); // durations of extraordinary events
				load_e[i] = gamma_p(generator2); // load of extraordinary events
			}
			else{
				T_e[i] = T_e[i-1] + exp_Te(generator2); // after an extra event, time until next extra event
				load_e[i] = 0;
			}
		}
		load_s[N] = gamma_s(generator2);
		if(load_e[N-1] == 0) load_e[N] = 0;
		else load_e[N] = gamma_p(generator2);
	}


	double operator() (double t){
		t = t / 8760; // express t in years
		double D_s = stepfun(t, T_s, load_s);
		double D_e = stepfun(t, T_e, load_e);
		double load_t = phi * R0 * (gamma * D_d + D_e + D_s)/(gamma * alpha_d + alpha_l);
		return load_t;
	}

	vector<double> max(){
		double max_load;
		vector<double> res(2,0);
		for(double t = t_start; t < t_end; t += constdt){
			max_load = (*this)(t);
			if(max_load > res[1]){
				res[0] = t;
				res[1] = max_load;
			}
		}
		return(res);
	}
/*
	void write(ofstream& file){
		if(file.is_open()){
			file << "T_s, load_s, T_e, load_e" << endl;
			for(int i = 0; i < T_s.size(); ++i){
				file << T_s[i] << "," << load_s[i] << "," << T_e[i] << "," << load_e[i] << endl;
			}
			file << " ," << load_s[T_s.size()] << ", ," << load_e[T_s.size()] << endl;
		}
		else{
			Rcout << "File is not open!" << endl;
		}
	}
*/
};

map<string, double> load_profile_par{
  {"phi", 1},
  {"load_s_shape", 3.122},
  {"load_s_scale", 0.0481},
  {"load_p_shape", 0.826},
  {"load_p_scale", 0.1023},
  {"R0", 2722},
  {"alpha_d", 1.25},
  {"alpha_l", 1.5},
  {"load_d_mean", 1},
  {"load_d_sd", 0.1},
  {"gamma", 0.25},
  {"mean_Ts", 10},
  {"mean_Te", 1},
  {"mean_Tp", 0.03835}
  //{"N", 100}
};

class constant_load{

	double load_level;

public:
	constant_load(double l): load_level(l) {};

	double operator() (double t){
		// t is in hour
		double T0 = load_level/k;
		if(t <= T0) return k*t;
		else return load_level;
	}

	vector<double> max(){
		double T0 = load_level/k;
		vector<double> res(2,0);
		res[0] = T0;
		res[1] = load_level;
		return(res);
	}

};

template <class T>
void CanADM(const double &x , double &dxdt , const double t, vector<double> parameter, double tau_s, T tau){

	double a, b ,c ,n, s0, mu;

	a = parameter[0];
	b = parameter[1];
	c = parameter[2];
	n = parameter[3];
	s0 = parameter[4];
	mu = parameter[6];

	dxdt = 1/mu * ( pow(a * tau_s * max((tau(t)/tau_s - s0), 0.0), b) + ( pow(c * tau_s * max((tau(t)/tau_s - s0), 0.0), n) * x));
}

double Teqn(double Ts, vector<double> params) {
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


class BadConversion : public std::runtime_error {
public:
	BadConversion(const std::string& s)
	: std::runtime_error(s)
	{ }
};

vector<double> generateC(vector<double> &theta){
	// Generate random variates
	normal_distribution<double> normA(theta[0], theta[1]);
	normal_distribution<double> normb(theta[2], theta[3]);
	normal_distribution<double> normC(theta[4], theta[5]);
	normal_distribution<double> normn(theta[6], theta[7]);
	normal_distribution<double> norms0(theta[8], theta[9]);

	vector<double> pars(7);
	pars[0] = exp(normA(generator2));
	pars[1] = exp(normb(generator2));
	pars[2] = exp(normC(generator2));
	pars[3] = exp(normn(generator2));
	double tempx = exp(norms0(generator2));
	pars[4] = tempx / (1.0+tempx);
	pars[5] = k;
	pars[6] = mu;

	return pars;
}


//[[Rcpp::export]]
void modifynPerTheta(const int newnPerTheta, const int verbose) {
  n_per_theta = newnPerTheta;
  if (verbose == 1) {
    Rcout << "number of simulated samples is set to be: " << n_per_theta << endl;
  }
}


//[[Rcpp::export]]
void modifyTEnd(const double newEndT, const int verbose) {
  t_end = newEndT;
  if (verbose == 1) {
    Rcout << "Return period is set to be: " << t_end << endl;
  }
}


//[[Rcpp::export]]
void modifyTimeStep(const double newdt, const int verbose) {
  constdt = newdt;
  if (verbose == 1) {
    Rcout << "Time step is set to be: " << constdt << endl;
  }
}

//' Display the current constants
//'
//' @description
//' Prints the constants currently associated with reliability assessment:
//'
//' * ramp load rate (k);
//'
//' * sample size of reliability simulation for each \eqn{\theta};
//'
//' * start time;
//'
//' * return period (in hours);
//'
//' * time step of ODE solver.
//'
//' @export
//'
//[[Rcpp::export]]
void catConstants() {
	//Rcout << "mu: " << mu << endl;
  Rcout << "standard ramp loading rate: " << k << endl;
  //Rcout << "load level: " << load_level << endl;
  //Rcout << "T0: " << T0 << endl;
	Rcout << "Number of simulated samples: " <<  n_per_theta << endl;
	Rcout << "Start period: " << t_start << endl;
  Rcout << "Return period: " << t_end << endl;
  Rcout << "Time step: " << constdt << endl;
}


//[[Rcpp::export]]
void modifyLoadP(std::string paramName, double param, const int verbose) {
  load_profile_par[paramName] = param;
  if (verbose == 1) {
    std::string res = paramName + " is set to be: ";
    Rcout << res << load_profile_par[paramName] << endl;
  }
}

//' Print the load profile parameters
//'
//' @description
//' Displays the value of parameters in the current load profile.
//'
//' @param paramName the name of the parameter in load profile
//'
//' @export
//'
//[[Rcpp::export]]
void print_LoadProfile(std::string paramName) {
  std::string res = paramName + " is ";
  Rcout << res << load_profile_par[paramName] << endl;
}


//[[Rcpp::export]]
List solveODE(NumericMatrix inputTheta, double phi, const int verbose) {
	// phi is set to be 1 by default
	double Ts;
	double tau_s;
	double x;
	double t;
	double x_max = 1;

	if (phi != 1) {
		load_profile_par["phi"] = phi;
		Rcout << "phi is set to be " << load_profile_par["phi"] << endl;
	}
	// TEST
	vector<vector<double>> theta;
	int nrow = inputTheta.rows();
	// copy each row of inputTheta into theta
	for (int i = 0; i < nrow; i++) {
		NumericVector curRow = inputTheta(i, _);
		vector<double> curVec = as<vector<double>>(curRow);
		theta.push_back(curVec);
	}

	int ntheta = theta.size();

	vector<int> nFail_per_theta(ntheta, 0);
	vector<int> nFail_per_theta_noDOL(ntheta, 0);
	vector<double> prop_Fail(ntheta, 0);
	vector<double> prop_Fail_noDOL(ntheta, 0);
	vector<double> time_to_failure(ntheta * n_per_theta);

	vector<double> pars;

	if (verbose == 1) {
	  Rcout << "Read " << ntheta << " thetas!" << endl;
	  Rcout << "Start simulation!" << endl;
	}

	// run simulation
  #ifdef _OPENMP
	  omp_set_num_threads(24);
  #endif
	#pragma omp parallel for schedule(dynamic) private(t, x, pars, Ts, tau_s)
	for(int i = 0; i < ntheta; ++i){
		int tmp = 0;
		int tmp_noDOL = 0;


		adams_bashforth_moulton< 5 , double > stepper;

		for(int j = 0; j < n_per_theta; ++j){



			pars = generateC(theta[i]);

			load_profile tau(load_profile_par);
			//constant_load tau(4500);

			try{
				pair<double, double> r = bisect([pars](double t) {return Teqn(t, pars);}, 0.001, 262980.0, tol);
				Ts = r.first + (r.second - r.first)/2;
			}
			catch(const boost::math::evaluation_error& e){
				time_to_failure[i * n_per_theta + j] = t_end/8760.0;
				continue;
			}

			tau_s = Ts * k;

			vector<double> max_t = tau.max();

			if(max_t[1] < tau_s * pars[4]){
				time_to_failure[i * n_per_theta + j] = t_end/8760.0;
			}
			else{

				// no DOL
				if(max_t[1] > tau_s){
					tmp_noDOL += 1;
					tmp += 1;
					continue;
				}

				// solve the ODE
				//CanADM model(pars, tau_s, tau);
				auto model = [&pars, &tau_s, &tau](const double &x , double &dxdt , const double t){
					CanADM(x , dxdt , t, pars, tau_s, tau);
				};


				x = 0; // initial state
				t = t_start;
				//times.push_back(t);
				//x_vec.push_back(x);
				stepper.initialize(model, x , t , constdt );
				t += constdt;
				while(t < t_end){
					stepper.do_step(model, x, t, constdt);
					//times.push_back(t);
					//x_vec.push_back(x);
					if( x > x_max || std::isnan(x)) // if x > 1, stop calculation
					{
						tmp += 1;
						break;
					}
					t += constdt;
				}

				time_to_failure[i * n_per_theta + j] = t/8760.0;


			}

			// the ode can be solved simply with the following; but we want it stops early when some criterion is met
			// so we have to solve it manually.

			//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
			//    model , x , t_start , t_end , dt, push_back_state_and_time(x_vec, times));

		}
		nFail_per_theta[i] = tmp;
		nFail_per_theta_noDOL[i] = tmp_noDOL;
		prop_Fail[i] = (double)nFail_per_theta[i] / (double) n_per_theta;
		prop_Fail_noDOL[i] = (double)nFail_per_theta_noDOL[i] / (double) n_per_theta;

		if (verbose == 1){
		  Rcout << i << " " << tmp << " " << tmp_noDOL << endl;
		}

	}
	// return
	NumericVector timeOut(time_to_failure.begin(), time_to_failure.end());
	//NEW
	NumericVector prop_Fail_Out(prop_Fail.begin(), prop_Fail.end());
	NumericVector prop_Fail_noDOL_Out(prop_Fail_noDOL.begin(), prop_Fail_noDOL.end());

	List result = List::create(timeOut,prop_Fail_Out,prop_Fail_noDOL_Out);
	return result;
}
