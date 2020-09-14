#include <cmath>
#include <Rcpp.h>
#include <math.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double dnorm_cpp(double x)
{
	return 1.0 / sqrt(2.0 * M_PI) * exp(-0.5 * x * x);
}


// [[Rcpp::export]]
double dnorm_d2_cpp(double x)
{
	return (x * x - 1.0) * dnorm_cpp(x);
}

// [[Rcpp::export]]
double fh_combt(double x, NumericVector t, NumericVector w, double h)
{
	const int k = static_cast<int>(t.size());
	double suma = 0;
	double dif = 0;
	for(int i = 0; i < k; ++i)
	{
		dif = (x - t[i]) / h;
		suma += w[i] * dnorm_cpp(dif);
	}
	suma /= h;

	return suma;
}

// [[Rcpp::export]]
double fh_combt_d2(double x, NumericVector t, NumericVector w, double h)
{
	const int k = static_cast<int>(t.size());
	double suma = 0;
	double dif = 0;
	for(int i = 0; i < k; ++i)
	{
		dif = (x - t[i]) / h;
		suma += w[i] * dnorm_d2_cpp(dif);
	}
	suma /= pow(h, 3.0);

	return suma;
}

// [[Rcpp::export]]
double curv_cpp(NumericVector t, NumericVector w, double h, double lim1, double lim2, int lgrid)
{
	const double cte = (lim2 - lim1) / lgrid;
	double fhtd2lim1 = fh_combt_d2(lim1, t, w, h);
	double fhtd2lim2 = fh_combt_d2(lim2, t, w, h);
	double suma = 0.5 * (fhtd2lim1 * fhtd2lim1 + fhtd2lim2 * fhtd2lim2);
	for(int i = 1; i < lgrid; ++i)
	{
		double xi = lim1 + i * cte;
		double fhtd2xi = fh_combt_d2(xi, t, w, h);
		suma += fhtd2xi * fhtd2xi;
	}
	suma *= cte;
	return suma;
}

// [[Rcpp::export]]
NumericVector dicoto_lambda(double lambda, int nith, double h0, double h1, double rho, NumericVector hist, NumericVector combt, NumericVector combw, double lim1, double lim2)
{

  NumericVector output(2), hseq(5), curv(5), objective(5);
  double newh0 = h0, newh1 = h1, newrho = rho;
  const int kcombt = static_cast<int>(combt.size());

  for(int it = 0; it < nith; ++it)
  {

    for(int i = 0; i < 5; ++i)
    {

      double auxhi = newh0 * pow(newrho, i);
      hseq[i] = auxhi;

      double suma = 0;
      double dif = 0;
      for(int ii = 0; ii < kcombt; ++ii)
      {
        double auxcombtii = combt[ii];
        dif = hist[ii] - fh_combt(auxcombtii, combt, combw, auxhi);
        suma += dif * dif;
      }

      double auxcurvi = curv_cpp(combt, combw, auxhi, lim1, lim2, 100);
      curv[i] = auxcurvi;
      objective[i] = suma + lambda * auxcurvi;

    }

    int minind = std::min_element(objective.begin(), objective.end()) - objective.begin();
    double gboot = hseq[minind];
    output[0] = gboot;
    output[1] = curv[minind];

    if(minind == 0)
    {
      newh1 = hseq[1];
      newh0 = gboot / newrho;
      newrho = exp((log(newh1)-log(newh0))/4);
    } else {
      if(minind == 4)
      {
        newh0 = hseq[3];
        newh1 = gboot * newrho;
        newrho = exp((log(newh1)-log(newh0))/4);
      } else {
        newh0 = hseq[minind - 1];
        newh1 = hseq[minind + 1];
        newrho = exp((log(newh1)-log(newh0))/4);
      }
    }

  }

  return output;

}



// [[Rcpp::export]]
double zeta_hist_p(NumericVector hist, NumericVector combt, NumericVector comby, NumericVector combw, double Af2_mixt, double l0, double l1, double h0, double h1, double lrho, double rho, int nitlambda, int nith, double lim1, double lim2)
{

  NumericVector lseq(5), dicoto_return(2), ldist(5), bw(5);
  double newl0 = l0, newl1 = l1, newlrho = lrho, gboot = 0;

  for(int itlambda = 0; itlambda < nitlambda; ++itlambda)
  {

    for(int i = 0; i < 5; ++i)
    {
      lseq[i] = newl0 * pow(newlrho, i);
    }

    for(int li = 0; li < 5; ++li)
    {
      double lambda = lseq[li];
      dicoto_return = dicoto_lambda(lambda, nith, h0, h1, rho, hist, combt, combw, lim1, lim2);
      bw[li] = dicoto_return[0];
      ldist[li] = std::abs(dicoto_return[1] - Af2_mixt);
    }

    int lminind = std::min_element(ldist.begin(), ldist.end()) - ldist.begin();
    gboot = bw[lminind];

    if(lminind == 0)
    {
      newl1 = lseq[1];
      newl0 = lseq[0] / newlrho;
      newlrho = exp((log(newl1)-log(newl0))/4);
    } else {
      if(lminind == 4)
      {
        newl0 = lseq[3];
        newl1 = lseq[4] * newlrho;
        newlrho = exp((log(newl1)-log(newl0))/4);
      } else {
        newl0 = lseq[lminind - 1];
        newl1 = lseq[lminind + 1];
        newlrho = exp((log(newl1)-log(newl0))/4);
      }
    }

  }

  return gboot;

}






// [[Rcpp::export]]
double pnorm_cpp_fun(double x)
{
  const double cte = 1.0/sqrt(2.0);
  return 0.5 * erfc(- x * cte);
}


// [[Rcpp::export]]
double dnorm_d1_cpp(double x)
{
	return -x * dnorm_cpp(x);
}

// [[Rcpp::export]]
double Fh_combt(double x, NumericVector t, NumericVector w, double h)
{
  const int k = static_cast<int>(t.size());
  double suma = 0;
  double dif = 0;
  for(int i = 0; i < k; ++i)
  {
    dif = (x - t[i]) / h;
    suma += w[i] * pnorm_cpp_fun(dif);
  }
  
  return suma;
}


// [[Rcpp::export]]
double fh_combt_d1(double x, NumericVector t, NumericVector w, double h)
{
	const int k = static_cast<int>(t.size());
	double suma = 0;
	double dif = 0;
	for(int i = 0; i < k; ++i)
	{
		dif = (x - t[i]) / h;
		suma += w[i] * dnorm_d1_cpp(dif);
	}
	suma /= h * h;

	return suma;
}



// [[Rcpp::export]]
double slope_cpp(NumericVector t, NumericVector w, double h, double lim1, double lim2, int lgrid)
{
	const double cte = (lim2 - lim1) / lgrid;
	double fhtd1lim1 = fh_combt_d1(lim1, t, w, h);
	double fhtd1lim2 = fh_combt_d1(lim2, t, w, h);
	double suma = 0.5 * (fhtd1lim1 * fhtd1lim1 + fhtd1lim2 * fhtd1lim2);
	for(int i = 1; i < lgrid; ++i)
	{
		double xi = lim1 + i * cte;
		double fhtd1xi = fh_combt_d1(xi, t, w, h);
		suma += fhtd1xi * fhtd1xi;
	}
	suma *= cte;
	return suma;
}


// [[Rcpp::export]]
NumericVector dicoto_lambda_dist(double lambda, int nith, double h0, double h1, double rho, NumericVector emp, NumericVector comby, NumericVector combt, NumericVector combw, double lim1, double lim2)
{

  NumericVector output(2), hseq(5), slope(5), objective(5);
  double newh0 = h0, newh1 = h1, newrho = rho;
  const int kcomby = static_cast<int>(comby.size());

  for(int it = 0; it < nith; ++it)
  {

    for(int i = 0; i < 5; ++i)
    {

      double auxhi = newh0 * pow(newrho, i);
      hseq[i] = auxhi;

      double suma = 0;
      double dif = 0;
      for(int ii = 0; ii < kcomby; ++ii)
      {
        double auxcombyii = comby[ii];
        dif = emp[ii] - Fh_combt(auxcombyii, combt, combw, auxhi);
        suma += dif * dif;
      }

      double auxslopei = slope_cpp(combt, combw, auxhi, lim1, lim2, 100);
      slope[i] = auxslopei;
      objective[i] = suma + lambda * auxslopei;

    }

    int minind = std::min_element(objective.begin(), objective.end()) - objective.begin();
    double gboot = hseq[minind];
    output[0] = gboot;
    output[1] = slope[minind];

    if(minind == 0)
    {
      newh1 = hseq[1];
      newh0 = gboot / newrho;
      newrho = exp((log(newh1)-log(newh0))/4);
    } else {
      if(minind == 4)
      {
        newh0 = hseq[3];
        newh1 = gboot * newrho;
        newrho = exp((log(newh1)-log(newh0))/4);
      } else {
        newh0 = hseq[minind - 1];
        newh1 = hseq[minind + 1];
        newrho = exp((log(newh1)-log(newh0))/4);
      }
    }

  }

  return output;

}




// [[Rcpp::export]]
double zeta_hist_p_dist(NumericVector emp, NumericVector combt, NumericVector comby, NumericVector combw, double Af1_mixt, double l0, double l1, double h0, double h1, double lrho, double rho, int nitlambda, int nith, double lim1, double lim2)
{

  NumericVector lseq(5), dicoto_return(2), ldist(5), bw(5);
  double newl0 = l0, newl1 = l1, newlrho = lrho, gboot = 0;

  for(int itlambda = 0; itlambda < nitlambda; ++itlambda)
  {

    for(int i = 0; i < 5; ++i)
    {
      lseq[i] = newl0 * pow(newlrho, i);
    }

    for(int li = 0; li < 5; ++li)
    {
      double lambda = lseq[li];
      dicoto_return = dicoto_lambda_dist(lambda, nith, h0, h1, rho, emp, comby, combt, combw, lim1, lim2);
      bw[li] = dicoto_return[0];
      ldist[li] = std::abs(dicoto_return[1] - Af1_mixt);
    }

    int lminind = std::min_element(ldist.begin(), ldist.end()) - ldist.begin();
    gboot = bw[lminind];

    if(lminind == 0)
    {
      newl1 = lseq[1];
      newl0 = lseq[0] / newlrho;
      newlrho = exp((log(newl1)-log(newl0))/4);
    } else {
      if(lminind == 4)
      {
        newl0 = lseq[3];
        newl1 = lseq[4] * newlrho;
        newlrho = exp((log(newl1)-log(newl0))/4);
      } else {
        newl0 = lseq[lminind - 1];
        newl1 = lseq[lminind + 1];
        newlrho = exp((log(newl1)-log(newl0))/4);
      }
    }

  }

  return gboot;

}
