#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;

void internal_buildLFactors(double * const LFactors, const int k, const NumericVector& t, const double h1, const double h2)
{
  const bool diff_h = (h1 != h2);

  NumericVector x(k);
  NumericVector * const x2 = diff_h ? new NumericVector(k) : NULL;

  for (int i = 0; i < k; i++) {

    // Computes (t[i]-t)/h argument to K2 and K4 functions
    for (int j = 0; j < k; j++) {
      const double tmp = (t[i] - t[j]);
      x[j] = tmp / h1;
      if (diff_h) {
        (*x2)[j] = tmp / h2;
      }
    }

    NumericVector d = dnorm(x);

    if (diff_h) {
      NumericVector d2 = dnorm(*x2);

      for (int j = 0; j < k; j++) {
        // (x^2-1)*dnorm(x)
        LFactors[i * k + j] = (x[j] * x[j] - 1.) * d[j];
        //  ((x^2-3)*(x^2-1)-2*x^2)*dnorm(x
        double tmp = ((*x2)[j]) * ((*x2)[j]);
        LFactors[i * k + j + k * k] = ((tmp - 3.) * (tmp - 1.) - 2. * tmp) * d2[j];
      }
    } else {
      for (int j = 0; j < k; j++) {
        const double x2 = x[j] * x[j]; // hides x2 ptr
        const double x2_1 = x2 - 1.;
        // (x^2-1)*dnorm(x)
        LFactors[i * k + j] = x2_1 * d[j];
        //  ((x^2-3)*(x^2-1)-2*x^2)*dnorm(x)
        LFactors[i * k + j + k * k] = ((x2 - 3.) * x2_1 - 2. * x2) * d[j];
      }
    }
  }

  if (diff_h) {
    delete x2;
  }
}

// [[Rcpp::export]]
NumericVector buildLFactors(NumericVector t, double h1, double h2)
{
  const int k = t.length();

  //[:,:,1] constains L1 and [:,:,2] contains L2
  NumericVector LFactors(Dimension(k, k, 2));

  internal_buildLFactors(&(LFactors[0]), k, t, h1, h2);

  return LFactors;
}


// [[Rcpp::export]]
NumericVector calcw_cpp(NumericVector xb, NumericVector y)
{ int r;

  const int n = xb.length();
  const int k = y.length() - 1;

  NumericVector wb(k);
  const double double_n = (double)n;

  // sapply(2:(leny-1),function(i)sum( x<y[i]&x>=y[i-1] ))
  for (int i = 0; i < (k - 1); i++) {
    const int iy = i + 1;
    r = 0;
    for (int j = 0; j < n; j++) {
      const double tmp = xb[j];
      r += (tmp < y[iy]) && (tmp >= y[iy - 1]);
    }
    wb[i] = r / double_n;
  }

  // sum(x<=y[leny]&x>=y[leny-1])
  r = 0;
  for (int j = 0; j < n; j++) {
    const double tmp = xb[j];
    r += (tmp <= y[k]) && (tmp >= y[k-1]);
  }
  wb[k-1] = r / double_n;

  return wb;
}

double internal_Window_helper(double * const xb,  const int n, const NumericVector& y,  const NumericVector& t,
                              const double * const Lfactors, double& J1_b, double& J2_b)
{ double mub = 0., sigmab = 0., J1tmp = 0., J2tmp = 0.;
  int r, xbi;

  // sort input xb
  std::sort(xb, xb + n);

  //const int n = xb.length();
  const int k = y.length() - 1;

  //double wb[k]; PROVOCA ERROR (CRAN)!!!
  double* wb = new double[k];
  const double double_n = (double)n;

  /*
   1) xb is sorted
   2) always y[i] <= y[i+1]
   3) xb can be < y[0] or > y[k+1]
   */

  // 1) skip xb < y[0]
  for (xbi = 0; (xbi < n) && (xb[xbi] < y[0]); xbi++);

  // 2) for i in [0..k-1], compute # of xbs / y[i] <= xb < y[i+1]
  for (int i = 0; i < (k - 1); i++) {
    r = xbi;
    while ((xbi < n) && (xb[xbi] < y[i+1])) {
      xbi++;
    }
    wb[i] = (xbi - r) / double_n;
  }

  r = xbi;
  while ((xbi < n) && (xb[xbi] <= y[k])) {
    xbi++;
  }
  wb[k - 1] = (xbi - r) / double_n;

  /*

  // sapply(2:(leny-1),function(i)sum( x<y[i]&x>=y[i-1] ))
  for (int i = 0; i < (k - 1); i++) {
    const int iy = i + 1;
    r = 0;
    for (int j = 0; j < n; j++) {
      const double tmp = xb[j];
      r += (tmp < y[iy]) && (tmp >= y[iy - 1]);
    }
    wb[i] = r / double_n;
  }

  // sum(x<=y[leny]&x>=y[leny-1])
  r = 0;
  for (int j = 0; j < n; j++) {
    const double tmp = xb[j];
    r += (tmp <= y[k]) && (tmp >= y[k-1]);
  }
  wb[k-1] = r / double_n;
  */

  // mub <- sum(wb*t)
  for (int i = 0; i < k; i++) {
    mub += wb[i] * t[i];
  }

  // sigmab[b] = sum(wb*(t-mub)^2)
  for (int i = 0; i < k; i++) {
    const double tmp = t[i] - mub;
    sigmab += wb[i] * (tmp * tmp);
  }

  //J1_b[b]   = L1factoredFunc(wb)
  //J2_b[b]   = L2factoredFunc(wb)
  for (int i = 0; i < k; i++) {
    double J1red = 0., J2red = 0.;
    for (int j = 0; j < k; j++) {
      J1red += Lfactors[i * k + j] * wb[j];
      J2red += Lfactors[i * k + j + k * k] * wb[j];
    }
    J1tmp += J1red * wb[i];
    J2tmp += J2red * wb[i];
  }
  J1_b = -J1tmp;
  J2_b = J2tmp;
  
  delete [] wb; // LINEA AÑADIDA (CRAN)

  return sqrt(sigmab);
}

// [[Rcpp::export]]
double Window_helper(int b, NumericVector xb, NumericVector y,
                     NumericVector t, NumericVector Lfactors,
                     NumericVector J1_b, NumericVector J2_b)
{
  return internal_Window_helper(&(xb[0]), xb.length(), y, t, &(Lfactors[0]), J1_b[b - 1], J2_b[b - 1]);
}


// [[Rcpp::export]]
void main_method_np(int hn, int B, NumericVector hseq, NumericMatrix xbm, NumericVector y,
                    NumericVector t, NumericVector MSE_J1, NumericVector MSE_J2, double J1_np, double J2_np, NumericVector J1b, NumericVector J2b)
{ double J1_b, J2_b;

//double J1mat[hn][B];
//double J2mat[hn][B];
double * const J1mat = new double[hn*B];
double * const J2mat = new double[hn*B];

  const int n = xbm.nrow();  // xbm is a (n x B) matrix
  const int k = t.length();

  double * const Lfactors = new double[k * k * 2];

  for (int hi = 0; hi < hn; hi++) {
    double avg_J1_b = 0.;
    double avg_J2_b = 0.;

    const double h = hseq[hi];
    const double h_pow_3 = h * h * h;
    const double h_pow_5 = h_pow_3 * h * h;

    //NumericVector Lfactors = buildLFactors(t, h, h);
    internal_buildLFactors(Lfactors, k, t, h, h);

    for (int b = 0; b < B; b++) {
      double *const xb = &(xbm[n*b]);
      const double sigmab = internal_Window_helper(xb, n, y, t, Lfactors, J1_b, J2_b);
      const double J1_b_f = (sigmab * sigmab * sigmab) * (J1_b / h_pow_3) - J1_np;
      const double J2_b_f = (sigmab * sigmab * sigmab * sigmab * sigmab) * (J2_b / h_pow_5) - J2_np;
	  //J1mat[hi][b] = J1_b_f+J1_np;
	  //J2mat[hi][b] = J2_b_f+J2_np;
	  J1mat[hi*B+b] = J1_b_f+J1_np;
	  J2mat[hi*B+b] = J2_b_f+J2_np;
      avg_J1_b += J1_b_f * J1_b_f;
      avg_J2_b += J2_b_f * J2_b_f;
    }
    MSE_J1[hi] = avg_J1_b / B;
    MSE_J2[hi] = avg_J2_b / B;
  }

  int i_h1 = 0;
  int i_h2 = 0;
  double mini1 = MSE_J1[0];
  double mini2 = MSE_J2[0];

  for(int i=0;i<hn;i++){
	if(MSE_J1[i] < mini1){
		mini1 = MSE_J1[i];
		i_h1 = i;
	}
	if(MSE_J2[i] < mini2){
		mini2 = MSE_J2[i];
		i_h1 = i;
	}
  }


  for(int i=0;i<B;i++){
	//J1b[i] = J1mat[i_h1][i];
	//J2b[i] = J2mat[i_h2][i];
	J1b[i] = J1mat[i_h1*B+i];
	J2b[i] = J2mat[i_h2*B+i];
  }

  delete [] Lfactors;
  delete [] J1mat;
  delete [] J2mat;
}
