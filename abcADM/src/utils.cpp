#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdlib.h> // for atof
#include <stdexcept> // for handling exceptions
#include "utils.h"


using namespace std;
const double PI  =3.141592653589793238463;


class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }
};

double convertToDouble(const string& s)
{
    istringstream i(s);
    double x;
    char c;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

vector< vector<double> > readCSV_to_vec(ifstream& file){
	vector< vector<double> > result;
	vector<double> row;
	double p;
	string line;
    string cell;

	while(getline(file, line)){

		stringstream lineStream(line);

		while(getline(lineStream, cell, ','))
		{
			p = convertToDouble(cell);
		    row.push_back(p);
		}
		result.push_back(row);
		row.clear();
	}
	return result;
}

vector<double> colMeans(vector< vector<double> > x){
	int ncol = x[1].size();
	int nrow = x.size();
	vector<double> res(ncol, 0);
	for(int i = 0; i < ncol; i++){
		for(int j = 0; j < nrow; j++){
			res[i] += x[j][i]/nrow;
		}
	}
	return res;
}

vector<double> colSds(vector< vector<double> > x){
	int ncol = x[1].size();
	int nrow = x.size();
	vector<double> res(ncol, 0);
	vector<double> colmean = colMeans(x);
	for(int i = 0; i < ncol; i++){
		for(int j = 0; j < nrow; j++){
			res[i] += pow(x[j][i], 2);
		}
		res[i] = (res[i] - pow(colmean[i], 2))/(nrow - 1);
	}
	return res;

}


double dist(vector<double> x, vector<double> y){
	double d = 0;
	int n = x.size();
	for(int i = 0; i < n; i++){
		d += pow((x[i] - y[i]), 2);
	}
	return sqrt(d);
}


vector<double> quantile(vector<double> x, vector<double> prob){
	sort(x.begin(), x.end());
	int n = x.size();
	vector<double> q;
	double h, temp;
	int lo, hi;
	for(auto p : prob){
		h = 1 + (n - 1) * p;
		lo = (int)floor(h);
		hi = (int)ceil(h);
		temp = x[lo-1] + (h - lo) * (x[hi-1] - x[lo-1]);
		q.push_back(temp);
	}
	return q;
}


double dinvgamma(double x, double shape, double scale, bool logp = false){
    if(x < 0) return(0);
    double log_density = shape * log(scale) - lgamma(shape) - (shape +
        1) * log(x) - (scale/x);
    if(logp) return(log_density);
    else return(exp(log_density));
}

double dnormC(double x, double mu, double sig, bool logp = false){
    double log_density = -log(sig * sqrt(2 * PI)) - pow((x - mu) / sig, 2) / 2;
    if(logp) return(log_density);
    else return(exp(log_density));
}


double log_prior(vector<double> theta, int n_par){
	double logp = 0;
	for(int i = 0; i < n_par-2; i+=2){
		logp += dnormC(theta[i], 0, 20, true);
	}
	logp += dnormC(theta[8], 0, 1, true); // stronger informative prior for sigma_0's mean (average 50% threshold)
	for(int i = 1; i < n_par; i+=2){
		logp += dinvgamma(pow(theta[i], 2), 0.01, 0.01, true);
	}
	return logp;
}
