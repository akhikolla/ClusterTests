#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "model.h"
#include "utils.h"

using namespace std;
using namespace Rcpp;

// global variables
const int n_par = 10; // number of parameters
const vector<double> probs = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}; // quantiles
const int n_s = probs.size(); // number of summary statistics
const vector<double> proposal_scale = {0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01, 0.1, 0.01};

// set up random number generator1
mt19937 generator1(time(0));

//[[Rcpp::export]]
List MCMC_ABC(const int numSamples, const int numBurning, const int numThining, const double delta,
	 const double k_val, const NumericVector wVec, const NumericVector censorVec, const List inputData) {
	/*
	Each element of inputData is a NumericVector
	*/
	const int m = numSamples; // number of ABC samples
	const int burning = numBurning;
	const int thining = numThining;
	const double k = k_val;

	int accept = 0;

	uniform_real_distribution<double> unif(0.0, 1.0);
	normal_distribution<double> norm_proposal(0.0, 1);

	vector<double> censored = as<vector<double>>(censorVec);
	vector<double> w0 = as<vector<double>>(wVec);


	////////////////////////////////////////////////////////////////////////////////
	int lstSize = inputData.size();
	vector<vector<double>> constdata(lstSize);
	for (int i = 0; i < lstSize; i++) {
		vector<double> vData = as<vector<double>>(inputData[i]);
		constdata[i] = vData;
	}
	//////////////////////////////////////////////////////////////////////////////// TEST constdata
	////////////////////////////  Test  ////////////////////////////////////////////
  /*
	int a = constdata.size();
	for(int i = 0; i < a; i++) {
		int b = constdata[i].size();
		cout << b << endl;
		for(int j = 0; j < b; j++) {
			cout << constdata[i][j] << " | ";
		}
		cout << endl;
	}
  */
////////////////////////////  Test  ////////////////////////////////////////////


	vector< vector<double> > s_obs(lstSize);
	vector<int> N(lstSize);
	vector<int> n_c_obs(lstSize);

	for(int i = 0; i < lstSize; ++i){
		s_obs[i] = summary_stat(constdata[i], censored[i], probs);
		N[i] = constdata[i].size();
		n_c_obs[i] = count_if(constdata[i].begin(), constdata[i].end(), [=](double c){return c > censored[i];});
	}

	vector<int> n_obs = N;


	vector< vector<double> > theta(m, vector<double>(n_par));
	vector< vector< vector<double> > > s(m, vector< vector<double> >(lstSize));
	vector< vector<double> > newData(lstSize);
	vector<double> newtheta(n_par), prevtheta(n_par);
	vector< vector<double> > s_simu(lstSize, vector<double>(n_s)), prev_s(lstSize, vector<double>(n_s));
	double u;
	double log_alpha;
	double log_den = 0;
	double log_num = 0;
	vector<int> n_censored(lstSize);



	// generate random samples
	theta[0] = {-8.33547,0.527351,3.77781,0.0567059,-24.2165,0.0897585,-0.847272,0.405729,-0.737646,0.179287}; // the initial value of theta

	log_den = log_prior(theta[0], n_par);
	for(int f = 0; f < lstSize; ++f){
		newData[f] = CANL(theta[0], N[f], w0[f], k, 0.5);
		s[0][f] = summary_stat(newData[f], censored[f], probs);
		n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [=](double c){return c > censored[f];});
		//log of the denominator in alpha
		log_den += dnormC(dist(s[0][f], s_obs[f])/delta, 0.0, 1.0, true);
		if (n_c_obs[f] > 0) {
			log_den += n_c_obs[f] * log((double)n_censored[f]/(double)N[f]) +
				(n_obs[f] - n_c_obs[f]) * log(1 - (double)n_censored[f]/(double)N[f]);
		}
	}



	// burning
	prevtheta = theta[0];
	prev_s = s[0];
	for(int i = 1; i < burning; ++i){

		// sample parameters from the proposal distribution
		for(int l = 0; l < n_par; ++l){
			newtheta[l] = prevtheta[l] + proposal_scale[l] * norm_proposal(generator1);
		}


		if(newtheta[1] < 0 || newtheta[3] < 0 || newtheta[5] < 0 || newtheta[7] < 0 || newtheta[9] < 0){
			continue;
		}

		// generate data
		log_num = log_prior(newtheta, n_par);
		for(int f = 0; f < lstSize; ++f){
			newData[f] = CANL(newtheta, N[f], w0[f], k, 0.5);
			s_simu[f] = summary_stat(newData[f], censored[f], probs);
			n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [=](double c){return c > censored[f];});
			//log of the numerator in alpha
			log_num += dnormC(dist(s_simu[f], s_obs[f])/delta, 0.0, 1.0, true);

			if (n_c_obs[f] > 0) {
					log_num+=	n_c_obs[f] * log((double)n_censored[f]/(double)N[f]) +
						(n_obs[f] - n_c_obs[f]) * log(1 - (double)n_censored[f]/(double)N[f]);
			}
			if (n_c_obs[f] == 0 && n_censored[f] > 0) {
			  log_num -= 10000.0;
			}
		}

		// accept/reject step
		u = unif(generator1);

		log_alpha = log_num - log_den;
		if(log(u) < log_alpha &&
			!std::isnan(log_alpha) &&
			all_of(newData.begin(), newData.end(), [](vector<double> d){return all_of(d.begin(), d.end(), [](double dd) {return dd > 0;});})){

			prevtheta = newtheta;
			log_den = log_num;
			prev_s = s_simu;
		}
	}


	theta[0] = prevtheta;
	s[0] = prev_s;

	for(int i = 1; i < m; ++i){

		for(int j = 0; j < thining; ++j){

			// sample parameters from the proposal distribution
			for(int l = 0; l < n_par; ++l){
				newtheta[l] = prevtheta[l] + proposal_scale[l] * norm_proposal(generator1);
			}

			if(newtheta[1] < 0 || newtheta[3] < 0 || newtheta[5] < 0 || newtheta[7] < 0 || newtheta[9] < 0){
				continue;
			}


			// generate data
			log_num = log_prior(newtheta, n_par);
			for(int f = 0; f < lstSize; ++f){
				newData[f] = CANL(newtheta, N[f], w0[f], k, 0.5);
				s_simu[f] = summary_stat(newData[f], censored[f], probs);
				n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [=](double c){return c > censored[f];});
				//log of the numerator in alpha
				log_num += dnormC(dist(s_simu[f], s_obs[f])/delta, 0.0, 1.0, true);
				if (n_c_obs[f] > 0) {
						log_num+=	n_c_obs[f] * log((double)n_censored[f]/(double)N[f]) +
							(n_obs[f] - n_c_obs[f]) * log(1 - (double)n_censored[f]/(double)N[f]);
				}
				if (n_c_obs[f] == 0 && n_censored[f] > 0) {
				  log_num -= 10000.0;
				}

			}

			// accept/reject step
			u = unif(generator1);

			log_alpha = log_num - log_den;
			if(log(u) < log_alpha &&
				!std::isnan(log_alpha) &&
				all_of(newData.begin(), newData.end(), [](vector<double> d){return all_of(d.begin(), d.end(), [](double dd) {return dd > 0;});})){

				prevtheta = newtheta;
				log_den = log_num;
				prev_s = s_simu;
			}

			if(log(u) < log_alpha &&
				!std::isnan(log_alpha) &&
				all_of(newData.begin(), newData.end(), [](vector<double> d){return all_of(d.begin(), d.end(), [](double dd) {return dd > 0;});})){

				prevtheta = newtheta;
				log_den = log_num;
				prev_s = s_simu;
				accept++;
			}
		}
		theta[i] = prevtheta;
		s[i] = prev_s;
	}

	//size of theta is n by n_par
	NumericMatrix outMat(m, n_par);
	for(int i = 0; i < m; i++) {
		NumericVector tmpOut(theta[i].begin(), theta[i].end());
		outMat(i, _) = tmpOut;
	}
	double acceptRate = (double)accept / (double)(m * thining);
	List result = List::create(acceptRate, outMat);
	return result;
}
