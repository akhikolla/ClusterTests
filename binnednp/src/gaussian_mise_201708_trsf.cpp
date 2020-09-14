#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif

namespace {

  struct MiseCapturedData {
    double d_n;
    int k;
    NumericVector *wg;
    NumericVector *w;
    NumericVector *t;
    double gboot;
    double AK;

    double *tdiff;

    MiseCapturedData() :
    wg(NULL),
    w(NULL),
    t(NULL),
    tdiff(NULL)
    {}
  };

  double sqrt2;

  MiseCapturedData Misecd;

  // https://github.com/SurajGupta/r-source/blob/master/src/nmath/dnorm.c
  // https://stackoverflow.com/questions/27425658/speed-up-dnorm-function
  // exp(-x^2/2)/sqrt(2*pi)
  // WARNING: ASSUMES mu and sigma are not NaNs
  void internal_dnorm_ms(double * const out, const double * const in, const int n, const double mu, const double sigma)
  { static const double dv = 1.0 / sqrt(2. * PI);

    //NumericVector nd_in(n);                //DBG
    //std::copy(in, in + n, nd_in.begin());  //DBG
    //NumericVector nd_out = dnorm(nd_in);   //DBG

    const double dv_sigma = dv / sigma;
    const double dv_num_sigma = 1.0 / ( 2. * sigma * sigma);

    for (int i = 0; i < n; ++i) {
      const double x = in[i] - mu;
      if (isnan(x)) {
        out[i] = x;
      } else {
        out[i] = exp(-(x * x) * dv_num_sigma) * dv_sigma;
      }
      //if (out[i] != nd_out[i]) {
      //  printf("%lf -> my %lf != %lf\n", x, out[i], nd_out[i]);
      //}
    }
  }

  /*
   void internal_dnorm(double * const out, const double * const in, const int n)
   {
   NumericVector nd_in(n);
   std::copy(in, in + n, nd_in.begin());
   NumericVector nd_out = dnorm(nd_in);
   std::copy(nd_out.begin(), nd_out.end(), out);
   }
   */

}

// [[Rcpp::export]]
void gaussian_mise_initialize(int n, int k, NumericVector wg, NumericVector w, NumericVector t, double gboot, double AK)
{

  if ( Misecd.wg != NULL) {
    delete  Misecd.wg;
    delete  Misecd.w;
    delete  Misecd.t;
    delete  [] Misecd.tdiff;
  } else {
    sqrt2 = sqrt(2.0); //only first time
  }

  Misecd.d_n = (double)n;
  Misecd.k = k;
  Misecd.wg = new NumericVector(wg);
  Misecd.w = new NumericVector(w);
  Misecd.t = new NumericVector(t);
  Misecd.gboot = gboot;
  Misecd.AK = AK;

  const int t_length = t.length(); // probably == k

  Misecd.tdiff = new double [k * t_length]; // k rows  x t_length columns

  for (int i = 0; i < k; i++) {
    double * const tmp = Misecd.tdiff + (t_length * i);
    const double t_i = t[i];
    for (int j = 0; j < t_length; j++) {
      tmp[j] = t_i - t[j];
    }
  }
}

// [[Rcpp::export]]
double gaussian_mise(double h)
{
  const double d_n = Misecd.d_n;
  const int k = Misecd.k;
  const NumericVector& wg = *Misecd.wg;
  const NumericVector& w  = *Misecd.w ;
  const NumericVector& t  = *Misecd.t ;
  const double gboot = Misecd.gboot;
  const double AK = Misecd.AK;

  const int t_length = t.length(); // probably == k

  //double d[t_length]; DA ERROR (CRAN)!!!
  double* d = new double[t_length];

  // suma1 <- sum( sapply( 1:k, function(i)sum( wg[i]*wg*dnorm(t[i]-t,0,h*sqrt(2)) ) ) )
  double suma1 = 0.;
  for (int i = 0; i < k; i++) {
    internal_dnorm_ms(d, Misecd.tdiff + (t_length * i), t_length, 0., h * sqrt2);
    double tmp = 0.;
    for (int j = 0; j < t_length; j++) {
      tmp += wg[j] * d[j];
    }
    suma1 += tmp * wg[i];
  }

  // suma2 <- sum( sapply( 1:k, function(i)sum( wg[i]*w*dnorm(t[i]-t,0,sqrt(h^2+gboot^2)) ) ) )
  double suma2 = 0.;
  for (int i = 0; i < k; i++) {
   internal_dnorm_ms(d, Misecd.tdiff + (t_length * i), t_length, 0., sqrt(h * h + gboot * gboot));
    double tmp = 0.;
    for (int j = 0; j < t_length; j++) {
      tmp += w[j] * d[j];
    }
    suma2 += tmp * wg[i];
  }

  // suma3 <- sum( sapply( 1:k, function(i)sum( w[i]*w*dnorm(t[i]-t,0,gboot*sqrt(2)) ) ) )
  double suma3 = 0.;
  for (int i = 0; i < k; i++) {
    internal_dnorm_ms(d, Misecd.tdiff + (t_length * i), t_length, 0., gboot * sqrt2);
    double tmp = 0.;
    for (int j = 0; j < t_length; j++) {
      tmp += w[j] * d[j];
    }
    suma3 += tmp * w[i];
  }
  
  delete [] d; // LINEA AÑADIDA (CRAN)!!!!

  return (d_n - 1.)/d_n * suma1 - 2. * suma2 + suma3 + AK / (d_n * h);
}

/*
// [[Rcpp::export]]
double gaussian_mise_loop(int hn, NumericVector seq, double rho)
{ double h2, seqnew[5];

  const double sqrt_rho = sqrt(rho);
  const int seq_length = 5; //seq.length();

  for (int i = 0; i < hn; ) {

    // j <- which.min(mise(seq))
    double min_val = gaussian_mise(seq[0]);
    int j = 0;
    for (int seq_index = 1; seq_index < seq_length; seq_index++) {
      const double tmp = gaussian_mise(seq[seq_index]);
      if (tmp < min_val) {
        min_val = tmp;
        j = seq_index;
      }
    }

    // h2 <- seq[j]
    h2 = seq[j];

    // i <- i + 1
    i++;

    const double tmp = pow(sqrt_rho, 1.0/i);
    switch (j) {
      case 0:
        // seq <- c(h2/rho^(2/i),h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),seq[2])
        seqnew[0] = h2 / pow(rho, 2.0/i);
        seqnew[1] = h2 / tmp;
        seqnew[2] = h2;
        seqnew[3] = h2 * tmp;
        seqnew[4] = seq[1];
        break;

      case 4:
        // seq <- c(seq[4],h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),h2*rho^(1/i))
        seqnew[0] = seq[3];
        seqnew[1] = h2 / tmp;
        seqnew[2] = h2;
        seqnew[3] = h2 * tmp;
        seqnew[4] = h2 * pow(rho, 1.0/i);
        break;

      default:

        // seq <- c(seq[j-1],h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),seq[j+1])
        seqnew[0] = seq[j-1];
        seqnew[1] = h2 / tmp;
        seqnew[2] = h2;
        seqnew[3] = h2 * tmp;
        seqnew[4] = seq[j+1];
        break;
    }

    for (int seq_index = 0; seq_index < seq_length; seq_index++) {
      seq[seq_index] = seqnew[seq_index];
    }
  }

  return h2;
}
*/

// [[Rcpp::export]]
double gaussian_mise_loop(int hn, NumericVector seq, double rho)
{ double h2, seqnew[5];

  //double sqrt_rho = sqrt(rho);
  const int seq_length = 5; //seq.length();

  for (int i = 0; i < hn; ) {

    // j <- which.min(mise(seq))
    double min_val = gaussian_mise(seq[0]);
    int j = 0;
    for (int seq_index = 1; seq_index < seq_length; seq_index++) {
      const double tmp = gaussian_mise(seq[seq_index]);
      if (tmp < min_val) {
        min_val = tmp;
        j = seq_index;
      }
    }

    // h2 <- seq[j]
    h2 = seq[j];

    // i <- i + 1
    i++;

    //const double tmp = pow(sqrt_rho, 1.0/i);
    switch (j) {
      case 0:
        // seq <- c(h2/rho^(2/i),h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),seq[2])
        seqnew[0] = h2 / rho;
        seqnew[4] = seq[1];
        rho = exp( (log(seqnew[4])-log(seqnew[0])) / 4. );
        seqnew[1] = seqnew[0] * rho;
        seqnew[2] = seqnew[0] * pow(rho,2);
        seqnew[3] = seqnew[0] * pow(rho,3);
        break;

      case 4:
        // seq <- c(seq[4],h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),h2*rho^(1/i))
        seqnew[0] = seq[3];
        seqnew[4] = h2 * rho;
        rho = exp( (log(seqnew[4])-log(seqnew[0])) / 4. );
        seqnew[1] = seqnew[0] * rho;
        seqnew[2] = seqnew[0] * pow(rho,2);
        seqnew[3] = seqnew[0] * pow(rho,3);
        break;

      default:

        // seq <- c(seq[j-1],h2/(sqrt(rho))^(1/i),h2,h2*sqrt(rho)^(1/i),seq[j+1])
        rho = exp( (log(seq[j+1]) - log(seq[j-1])) / 4. );
        seqnew[0] = seq[j-1];
        seqnew[1] = seqnew[0] * rho;
        seqnew[2] = seqnew[0] * pow(rho,2);
        seqnew[3] = seqnew[0] * pow(rho,3);
        seqnew[4] = seq[j+1];
        break;
    }

    for (int seq_index = 0; seq_index < seq_length; seq_index++) {
      seq[seq_index] = seqnew[seq_index];
    }
  }

  return h2;
}

// [[Rcpp::export]]
double gaussian_dichotomy(int hn, NumericVector t)
{
  const int t_length = t.length();

  double min_t = t[0];
  double max_t = t[0];
  double min_diff = t[1] - t[0];

  for (int i = 1; i < t_length; i++) {
    const double tmp = t[i];
    if(tmp < min_t) min_t = tmp;
    if (tmp > max_t) max_t = tmp;
    const double diff = tmp - t[i-1];
    if(((diff < min_diff) && (diff > 1e-6)) || (min_diff < 1e-6)) min_diff = diff;
  }

  const double h0 = min_diff / 2.;
  const double h1 = max_t - min_t;
  //double h0 = min_diff / 2.;
  //double h1 = max_t - min_t;
  const int r = 5;
  //const double rho = exp( (log(h1)-log(h0)) / (r-1) );
  double rho = exp( (log(h1)-log(h0)) / (r-1) );

  NumericVector seq(r);
  for (int seq_index = 0; seq_index < r; seq_index++) {
    seq[seq_index] = h0 * pow(rho, (double)seq_index);
  }

  return gaussian_mise_loop(hn, seq, rho);
}
