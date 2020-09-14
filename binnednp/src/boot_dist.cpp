#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double pnorm_cpp(double x)
{
	const double cte = 1.0/sqrt(2.0);
	return 0.5 * erfc(- x * cte);
}


// [[Rcpp::export]]
double Fg(double x, NumericVector w, NumericVector t, double g)
{
	const int k = static_cast<int>(t.size());
	double suma = 0;
	for(int i = 0; i < k; i++)
	{
		double arg = (x - t[i]) / g;
		suma += w[i] * pnorm_cpp(arg);
	}
	return suma;
}


// [[Rcpp::export]]
double biasFh(double x, int n, NumericVector t, NumericVector w, NumericVector p, double g, double h)
{

	const int k = static_cast<int>(t.size());
	double suma = 0;
	for(int i = 0; i < k; i++)
	{
		double arg = (x - t[i]) / h;
		suma += pnorm_cpp(arg) * p[i];
	}
	suma += -Fg(x, w, t, g);
	return suma;

}


// [[Rcpp::export]]
double varFh(double x, int n, NumericVector t, NumericVector p, double h)
{

	const int k = static_cast<int>(t.size());
	double invn = 1.0/n;
	double suma1 = 0;
	double suma2 = 0;
	for(int i = 0; i < k; i++)
	{
		double arg = (x - t[i]) / h;
		double pnormarg = pnorm_cpp(arg);
		suma1 += pnormarg * pnormarg * p[i] * (1 - p[i]);
		for(int j = i+1; j < k; j++)
		{
			double arg2 = (x - t[j]) / h;
			double pnormarg2 = pnorm_cpp(arg2);
			suma2 += pnormarg * pnormarg2 * p[i] * p[j];
		}
	}
	return invn * suma1 - 2.0 * invn * suma2;

}



// [[Rcpp::export]]
double mise_Fh(double h, int n, NumericVector t, NumericVector w, NumericVector p, double g, int lgrid, double lim1, double lim2)
{
	const double cte = (lim2 - lim1) / lgrid;
	double biasFhlim1 = biasFh(lim1, n, t, w, p, g, h);
	double biasFhlim2 = biasFh(lim2, n, t, w, p, g, h);
	double varFhlim1 = varFh(lim1, n, t, p, h);
	double varFhlim2 = varFh(lim2, n, t, p, h);
	double suma = 0.5 * (biasFhlim1 * biasFhlim1 + varFhlim1 + biasFhlim2 * biasFhlim2 + varFhlim2);
	for(int i = 1; i < lgrid-1; i++)
	{
		double xi = lim1 + i * cte;
		double biasFhxi = biasFh(xi, n, t, w, p, g, h);
		double varFhxi = varFh(xi, n, t, p, h);
		suma += biasFhxi * biasFhxi + varFhxi;
	}
	suma *= cte;
	return suma;
}


// [[Rcpp::export]]
double boot_bw_dist(int nit, double h0, double h1, double rho, int n, NumericVector t, NumericVector w, NumericVector p, double g, int lgrid, double lim1, double lim2)
{
	NumericVector mises(5), hseq(5);
	int j = 2;
	double newh0 = h0, newh1 = h1, newrho = rho;

	for(int it = 0; it < nit; it++)
	{

		for(int i = 0; i < 5; i++)
		{
			hseq[i] = newh0 * pow(newrho, i);
			mises[i] = mise_Fh(hseq[i], n, t, w, p, g, lgrid, lim1, lim2);
		}
		j = std::min_element(mises.begin(), mises.end()) - mises.begin();

		if(j == 0) {
			newh0 = hseq[0] / newrho;
			newh1 = hseq[1];
			newrho = exp((log(newh1) - log(newh0)) / 4);
		} else {
			if(j == 4) {
				newh0 = hseq[3];
				newh1 = hseq[4] * newrho;
				exp((log(newh1) - log(newh0)) / 4);
			} else {
				newh0 = hseq[j - 1];
				newh1 = hseq[j + 1];
				newrho = exp((log(newh1) - log(newh0)) / 4);
			}
		}

	}

	return hseq[j];
}
