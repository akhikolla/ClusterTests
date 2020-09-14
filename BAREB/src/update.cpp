#include "update.h"

//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
  //
  //  Discussion:
  //
  //    The pseudorandom number should be uniformly distributed
  //    between A and B.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    02 October 2012
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Paul Bratley, Bennett Fox, Linus Schrage,
  //    A Guide to Simulation,
  //    Second Edition,
  //    Springer, 1987,
  //    ISBN: 0387964673,
  //    LC: QA76.9.C65.B73.
  //
  //    Bennett Fox,
  //    Algorithm 647:
  //    Implementation and Relative Efficiency of Quasirandom
  //    Sequence Generators,
  //    ACM Transactions on Mathematical Software,
  //    Volume 12, Number 4, December 1986, pages 362-376.
  //
  //    Pierre L'Ecuyer,
  //    Random Number Generation,
  //    in Handbook of Simulation,
  //    edited by Jerry Banks,
  //    Wiley, 1998,
  //    ISBN: 0471134031,
  //    LC: T57.62.H37.
  //
  //    Peter Lewis, Allen Goodman, James Miller,
  //    A Pseudo-Random Number Generator for the System/360,
  //    IBM Systems Journal,
  //    Volume 8, Number 2, 1969, pages 136-143.
  //
  //  Parameters:
  //
  //    Input, int A, B, the limits of the interval.
  //
  //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
  //    On output, SEED has been updated.
  //
  //    Output, int I4_UNIFORM_AB, a number between A and B.
  //
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;


  //
  //  Guarantee A <= B.
  //
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
  //
  //  Scale R to lie between A-0.5 and B+0.5.
  //
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 )
    +         r   * ( ( float ) b + 0.5 );
  //
  //  Use rounding to convert R to an integer between A and B.
  //
  value = round ( r );
  //
  //  Guarantee A <= VALUE <= B.
  //
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

double normal_01_cdf ( double x )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_CDF evaluates the Normal 01 CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    10 February 1999
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    A G Adams,
  //    Areas Under the Normal Curve,
  //    Algorithm 39,
  //    Computer j.,
  //    Volume 12, pages 197-198, 1969.
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the CDF.
  //
  //    Output, double CDF, the value of the CDF.
  //
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
  //
  //  |X| <= 1.28.
  //
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
                                                               + a6 / ( y + a7 ) ) ) );
    //
    //  1.28 < |X| <= 12.7
    //
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1
                               + b2  / ( fabs ( x ) + b3
                                           + b4  / ( fabs ( x ) - b5
                                           + b6  / ( fabs ( x ) + b7
                                           - b8  / ( fabs ( x ) + b9
                                           + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
                                           //
                                           //  12.7 < |X|
                                           //
  }
  else
  {
    q = 0.0;
  }
  //
  //  Take account of negative X.
  //
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}
//****************************************************************************80

double normal_01_cdf_inv ( double p )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_CDF_INV inverts the standard normal CDF.
  //
  //  Discussion:
  //
  //    The result is accurate to about 1 part in 10**16.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    27 December 2004
  //
  //  Author:
  //
  //    Original FORTRAN77 version by Michael Wichura.
  //    C++ version by John Burkardt.
  //
  //  Reference:
  //
  //    Michael Wichura,
  //    The Percentage Points of the Normal Distribution,
  //    Algorithm AS 241,
  //    Applied Statistics,
  //    Volume 37, Number 3, pages 477-484, 1988.
  //
  //  Parameters:
  //
  //    Input, double P, the value of the cumulative probability
  //    densitity function.  0 < P < 1.  If P is outside this range, an
  //    "infinite" value is returned.
  //
  //    Output, double NORMAL_01_CDF_INV, the normal deviate value
  //    with the property that the probability of a standard normal deviate being
  //    less than or equal to this value is P.
  //
{
  double a[8] = {
    3.3871328727963666080,     1.3314166789178437745E+2,
    1.9715909503065514427E+3,  1.3731693765509461125E+4,
    4.5921953931549871457E+4,  6.7265770927008700853E+4,
    3.3430575583588128105E+4,  2.5090809287301226727E+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252E+1,
    6.8718700749205790830E+2,  5.3941960214247511077E+3,
    2.1213794301586595867E+4,  3.9307895800092710610E+4,
    2.8729085735721942674E+4,  5.2264952788528545610E+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770E-1,
    2.27238449892691845833E-2,  7.74545014278341407640E-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550E-1,
    1.48103976427480074590E-1,  1.51986665636164571966E-2,
    5.47593808499534494600E-4,  1.05075007164441684324E-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230E-1,
    2.65321895265761230930E-2,  1.24266094738807843860E-3,
    2.71155556874348757815E-5,  2.01033439929228813265E-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690E-1,
    1.36929880922735805310E-1,  1.48753612908506148525E-2,
    7.86869131145613259100E-4,  1.84631831751005468180E-5,
    1.42151175831644588870E-7,  2.04426310338993978564E-15 };
  double q;
  double r;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;

  if ( p <= 0.0 )
  {
    value = -r8_huge ( );
    return value;
  }

  if ( 1.0 <= p )
  {
    value = r8_huge ( );
    return value;
  }

  q = p - 0.5;

  if ( fabs ( q ) <= split1 )
  {
    r = const1 - q * q;
    value = q * r8poly_value_horner ( 7, a, r )
      / r8poly_value_horner ( 7, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = r8_huge ( );
    }
    else
    {
      r = sqrt ( - log ( r ) );

      if ( r <= split2 )
      {
        r = r - const2;
        value = r8poly_value_horner ( 7, c, r )
          / r8poly_value_horner ( 7, d, r );
      }
      else
      {
        r = r - split2;
        value = r8poly_value_horner ( 7, e, r )
          / r8poly_value_horner ( 7, f, r );
      }
    }

    if ( q < 0.0 )
    {
      value = - value;
    }

  }

  return value;
}
//****************************************************************************80

void normal_01_cdf_values ( int &n_data, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
  //
  //  Discussion:
  //
  //    In Mathematica, the function can be evaluated by:
  //
  //      Needs["Statistics`ContinuousDistributions`"]
  //      dist = NormalDistribution [ 0, 1 ]
  //      CDF [ dist, x ]
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    28 August 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 17

  static double fx_vec[N_MAX] = {
    0.5000000000000000E+00,
    0.5398278372770290E+00,
    0.5792597094391030E+00,
    0.6179114221889526E+00,
    0.6554217416103242E+00,
    0.6914624612740131E+00,
    0.7257468822499270E+00,
    0.7580363477769270E+00,
    0.7881446014166033E+00,
    0.8159398746532405E+00,
    0.8413447460685429E+00,
    0.9331927987311419E+00,
    0.9772498680518208E+00,
    0.9937903346742239E+00,
    0.9986501019683699E+00,
    0.9997673709209645E+00,
    0.9999683287581669E+00 };

  static double x_vec[N_MAX] = {
    0.0000000000000000E+00,
    0.1000000000000000E+00,
    0.2000000000000000E+00,
    0.3000000000000000E+00,
    0.4000000000000000E+00,
    0.5000000000000000E+00,
    0.6000000000000000E+00,
    0.7000000000000000E+00,
    0.8000000000000000E+00,
    0.9000000000000000E+00,
    0.1000000000000000E+01,
    0.1500000000000000E+01,
    0.2000000000000000E+01,
    0.2500000000000000E+01,
    0.3000000000000000E+01,
    0.3500000000000000E+01,
    0.4000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double normal_01_mean ( )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    18 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Output, double MEAN, the mean of the PDF.
  //
{
  double mean;

  mean = 0.0;

  return mean;
}
//****************************************************************************80

double normal_01_moment ( int order )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_MOMENT evaluates moments of the Normal 01 PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    31 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Output, double NORMAL_01_MOMENT, the value of the moment.
  //
{
  double value;

  if ( ( order % 2 ) == 0 )
  {
    value = r8_factorial2 ( order - 1 );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double normal_01_pdf ( double x )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_PDF evaluates the Normal 01 PDF.
  //
  //  Discussion:
  //
  //    The Normal 01 PDF is also called the "Standard Normal" PDF, or
  //    the Normal PDF with 0 mean and standard deviation 1.
  //
  //    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    18 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the PDF.
  //
  //    Output, double PDF, the value of the PDF.
  //
{
  double pdf;
  const double r8_pi = 3.14159265358979323;

  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * r8_pi );

  return pdf;
}
//****************************************************************************80

double normal_01_sample ( int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_SAMPLE samples the standard normal probability distribution.
  //
  //  Discussion:
  //
  //    The standard normal probability distribution function (PDF) has
  //    mean 0 and standard deviation 1.
  //
  //    The Box-Muller method is used, which is efficient, but
  //    generates two values at a time.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    26 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input/output, int &SEED, a seed for the random number generator.
  //
  //    Output, double NORMAL_01_SAMPLE, a normally distributed random value.
  //
{
  double r1;
  double r2;
  const double r8_pi = 3.14159265358979323;
  double x;

  r1 = r8_uniform_01 ( seed );
  r2 = r8_uniform_01 ( seed );
  x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}
//****************************************************************************80

double normal_01_variance ( )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    18 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Output, double VARIANCE, the variance of the PDF.
  //
{
  double variance;

  variance = 1.0;

  return variance;
}
//****************************************************************************80

double normal_ms_cdf ( double x, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_CDF evaluates the Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    19 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the CDF.
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Output, double NORMAL_MS_CDF, the value of the CDF.
  //
{
  double cdf;
  double y;

  y = ( x - mu ) / sigma;

  cdf = normal_01_cdf ( y );

  return cdf;
}
//****************************************************************************80

double normal_ms_cdf_inv ( double cdf, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_CDF_INV inverts the Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    19 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Algorithm AS 111,
  //    Applied Statistics,
  //    Volume 26, pages 118-121, 1977.
  //
  //  Parameters:
  //
  //    Input, double CDF, the value of the CDF.
  //    0.0 <= CDF <= 1.0.
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Output, double NORMAL_MS_CDF_INV, the corresponding argument.
  //
{
  double x;
  double x2;

  x2 = normal_01_cdf_inv ( cdf );

  x = mu + sigma * x2;

  return x;
}
//****************************************************************************80

double normal_ms_mean ( double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_MEAN returns the mean of the Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    19 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Output, double NORMAL_MS_MEAN, the mean of the PDF.
  //
{
  double mean;

  mean = mu;

  return mean;
}
//****************************************************************************80

double normal_ms_moment ( int order, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_MOMENT evaluates moments of the Normal PDF.
  //
  //  Discussion:
  //
  //    The formula was posted by John D Cook.
  //
  //    Order  Moment
  //    -----  ------
  //      0    1
  //      1    mu
  //      2    mu^2 +         sigma^2
  //      3    mu^3 +  3 mu   sigma^2
  //      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
  //      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
  //      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
  //      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
  //      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6
  //           + 105 sigma^8
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    31 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Input, double MU, the mean of the distribution.
  //
  //    Input, double SIGMA, the standard deviation of the distribution.
  //
  //    Output, double NORMAL_MS_MOMENT, the value of the central moment.
  //
{
  int j;
  int j_hi;
  double value;

  j_hi = ( order / 2 );

  value = 0.0;
  for ( j = 0; j <= j_hi; j++ )
  {
    value = value
    + r8_choose ( order, 2 * j )
    * r8_factorial2 ( 2 * j - 1 )
    * pow ( mu, order - 2 * j ) * pow ( sigma, 2 * j );
  }

  return value;
}
//****************************************************************************80

double normal_ms_moment_central ( int order, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_MOMENT_CENTRAL evaluates central moments of the Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    31 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Input, double MU, the mean of the distribution.
  //
  //    Input, double SIGMA, the standard deviation of the distribution.
  //
  //    Output, double NORMAL_MS_MOMENT_CENTRAL, the value of the central moment.
  //
{
  double value;

  if ( ( order % 2 ) == 0 )
  {
    value = r8_factorial2 ( order - 1 ) * pow ( sigma, order );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double normal_ms_moment_central_values ( int order, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_MOMENT_CENTRAL_VALUES: moments 0 through 10 of the Normal PDF.
  //
  //  Discussion:
  //
  //    The formula was posted by John D Cook.
  //
  //    Order  Moment
  //    -----  ------
  //      0    1
  //      1    0
  //      2    sigma^2
  //      3    0
  //      4    3 sigma^4
  //      5    0
  //      6    15 sigma^6
  //      7    0
  //      8    105 sigma^8
  //      9    0
  //     10    945 sigma^10
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    31 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER <= 10.
  //
  //    Input, double MU, the mean of the distribution.
  //
  //    Input, double SIGMA, the standard deviation of the distribution.
  //
  //    Output, double NORMAL_MS_MOMENT_CENTRAL_VALUES, the value of the
  //    central moment.
  //
{
  double value;

  if ( order == 0 )
  {
    value = 1.0;
  }
  else if ( order == 1 )
  {
    value = 0.0;
  }
  else if ( order == 2 )
  {
    value = pow ( sigma, 2 );
  }
  else if ( order == 3 )
  {
    value = 0.0;
  }
  else if ( order == 4 )
  {
    value = 3.0 * pow ( sigma, 4 );
  }
  else if ( order == 5 )
  {
    value = 0.0;
  }
  else if ( order == 6 )
  {
    value = 15.0 * pow ( sigma, 6 );
  }
  else if ( order == 7 )
  {
    value = 0.0;
  }
  else if ( order == 8 )
  {
    value = 105.0 * pow ( sigma, 8 );
  }
  else if ( order == 9 )
  {
    value = 0.0;
  }
  else if ( order == 10 )
  {
    value = 945.0 * pow ( sigma, 10 );
  }


  return value;
}
//****************************************************************************80

double normal_ms_moment_values ( int order, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_MOMENT_VALUES evaluates moments 0 through 8 of the Normal PDF.
  //
  //  Discussion:
  //
  //    The formula was posted by John D Cook.
  //
  //    Order  Moment
  //    -----  ------
  //      0    1
  //      1    mu
  //      2    mu^2 +         sigma^2
  //      3    mu^3 +  3 mu   sigma^2
  //      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
  //      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
  //      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
  //      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
  //      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6
  //           + 105 sigma^8
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    31 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER <= 8.
  //
  //    Input, double MU, the mean of the distribution.
  //
  //    Input, double SIGMA, the standard deviation of the distribution.
  //
  //    Output, double NORMAL_MS_MOMENT_VALUES, the value of the central moment.
  //
{
  double value;

  if ( order == 0 )
  {
    value = 1.0;
  }
  else if ( order == 1 )
  {
    value = mu;
  }
  else if ( order == 2 )
  {
    value = pow ( mu, 2 ) + pow ( sigma, 2 );
  }
  else if ( order == 3 )
  {
    value = pow ( mu, 3 ) + 3.0 * mu * pow ( sigma, 2 );
  }
  else if ( order == 4 )
  {
    value = pow ( mu, 4 ) + 6.0 * pow ( mu, 2 ) * pow ( sigma, 2 )
    + 3.0 * pow ( sigma, 4 );
  }
  else if ( order == 5 )
  {
    value = pow ( mu, 5 ) + 10.0 * pow ( mu, 3 ) * pow ( sigma, 2 )
    + 15.0 * mu * pow ( sigma, 4 );
  }
  else if ( order == 6 )
  {
    value = pow ( mu, 6 ) + 15.0 * pow ( mu, 4 ) * pow ( sigma, 2 )
    + 45.0 * pow ( mu, 2 ) * pow ( sigma, 4 )
    + 15.0 * pow ( sigma, 6 );
  }
  else if ( order == 7 )
  {
    value = pow ( mu, 7 ) + 21.0 * pow ( mu, 5 ) * pow ( sigma, 2 )
    + 105.0 * pow ( mu, 3 ) * pow ( sigma, 4 )
    + 105.0 * mu * pow ( sigma, 6 );
  }
  else if ( order == 8 )
  {
    value = pow ( mu, 8 ) + 28.0 * pow ( mu, 6 ) * pow ( sigma, 2 )
    + 210.0 * pow ( mu, 4 ) * pow ( sigma, 4 )
    + 420.0 * pow ( mu, 2 ) * pow ( sigma, 6 ) + 105.0 * pow ( sigma, 8 );
  }


  return value;
}
//****************************************************************************80

double normal_ms_pdf ( double x, double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_PDF evaluates the Normal PDF.
  //
  //  Discussion:
  //
  //    PDF(MU,SIGMA;X)
  //      = exp ( - 0.5 * ( ( X - MU ) / SIGMA )^2 )
  //      / ( SIGMA * SQRT ( 2 * PI ) )
  //
  //    The normal PDF is also known as the Gaussian PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    19 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the PDF.
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Output, double NORMAL_MS_PDF, the value of the PDF.
  //
{
  double pdf;
  const double r8_pi = 3.14159265358979323;
  double y;

  y = ( x - mu ) / sigma;

  pdf = exp ( - 0.5 * y * y )  / ( sigma * sqrt ( 2.0 * r8_pi ) );

  return pdf;
}
//****************************************************************************80

double normal_ms_sample ( double mu, double sigma, int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_SAMPLE samples the Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    20 November 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Input/output, int &SEED, a seed for the random number generator.
  //
  //    Output, double NORMAL_MS_SAMPLE, a sample of the PDF.
  //
{
  double x;

  x = normal_01_sample ( seed );

  x = mu + sigma * x;

  return x;
}
//****************************************************************************80

double normal_ms_variance ( double mu, double sigma )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    NORMAL_MS_VARIANCE returns the variance of the Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    19 September 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the parameters of the PDF.
  //    0.0 < SIGMA.
  //
  //    Output, double NORMAL_MS_VARIANCE, the variance of the PDF.
  //
{
  double variance;

  variance = sigma * sigma;

  return variance;
}
//****************************************************************************80

double r8_abs ( double x )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_ABS returns the absolute value of an R8.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 November 2006
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the quantity whose absolute value is desired.
  //
  //    Output, double R8_ABS, the absolute value of X.
  //
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_choose ( int n, int k )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
  //
  //  Discussion:
  //
  //    The value is calculated in such a way as to avoid overflow and
  //    roundoff.  The calculation is done in R8 arithmetic.
  //
  //    The formula used is:
  //
  //      C(N,K) = N! / ( K! * (N-K)! )
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    09 June 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    ML Wolfson, HV Wright,
  //    Algorithm 160:
  //    Combinatorial of M Things Taken N at a Time,
  //    Communications of the ACM,
  //    Volume 6, Number 4, April 1963, page 161.
  //
  //  Parameters:
  //
  //    Input, int N, K, the values of N and K.
  //
  //    Output, double R8_CHOOSE, the number of combinations of N
  //    things taken K at a time.
  //
{
  int i;
  int mn;
  int mx;
  double value;

  if ( k < n - k )
  {
    mn = k;
    mx = n - k;
  }
  else
  {
    mn = n - k;
    mx = k;
  }

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_factorial2 ( int n )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_FACTORIAL2 computes the double factorial function.
  //
  //  Discussion:
  //
  //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
  //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
  //
  //  Example:
  //
  //     N Value
  //
  //     0     1
  //     1     1
  //     2     2
  //     3     3
  //     4     8
  //     5    15
  //     6    48
  //     7   105
  //     8   384
  //     9   945
  //    10  3840
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    22 January 2008
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the argument of the double factorial
  //    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
  //
  //    Output, double R8_FACTORIAL2, the value of Factorial2(N).
  //
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
//****************************************************************************80

void r8_factorial2_values ( int &n_data, int &n, double &f )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_FACTORIAL2_VALUES returns values of the double factorial function.
  //
  //  Formula:
  //
  //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
  //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
  //
  //    In Mathematica, the function can be evaluated by:
  //
  //      n!!
  //
  //  Example:
  //
  //     N    N!!
  //
  //     0     1
  //     1     1
  //     2     2
  //     3     3
  //     4     8
  //     5    15
  //     6    48
  //     7   105
  //     8   384
  //     9   945
  //    10  3840
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    07 February 2015
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //    Daniel Zwillinger,
  //    CRC Standard Mathematical Tables and Formulae,
  //    30th Edition,
  //    CRC Press, 1996, page 16.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, int &N, the argument of the function.
  //
  //    Output, double &F, the value of the function.
  //
{
# define N_MAX 16

  static double f_vec[N_MAX] = {
    1.0,
    1.0,
    2.0,
    3.0,
    8.0,
    15.0,
    48.0,
    105.0,
    384.0,
    945.0,
    3840.0,
    10395.0,
    46080.0,
    135135.0,
    645120.0,
    2027025.0 };

  static int n_vec[N_MAX] = {
    0,
    1,  2,  3,  4,  5,
    6,  7,  8,  9, 10,
    11, 12, 13, 14, 15 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    f = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double r8_huge ( )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_HUGE returns a "huge" R8.
  //
  //  Discussion:
  //
  //    HUGE_VAL is the largest representable legal real number, and is usually
  //    defined in math.h, or sometimes in stdlib.h.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    08 May 2003
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Output, double R8_HUGE, a "huge" real value.
  //
{
  return HUGE_VAL;
}
//****************************************************************************80

double r8_log_2 ( double x )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_LOG_2 returns the logarithm base 2 of the absolute value of an R8.
  //
  //  Discussion:
  //
  //    value = Log ( |X| ) / Log ( 2.0 )
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    22 March 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the number whose base 2 logarithm is desired.
  //    X should not be 0.
  //
  //    Output, double R8_LOG_2, the logarithm base 2 of the absolute
  //    value of X.  It should be true that |X| = 2^R_LOG_2.
  //
{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( 2.0 );
  }

  return value;
}
//****************************************************************************80

double r8_mop ( int i )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_MOP returns the I-th power of -1 as an R8 value.
  //
  //  Discussion:
  //
  //    An R8 is an double value.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    16 November 2007
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int I, the power of -1.
  //
  //    Output, double R8_MOP, the I-th power of -1.
  //
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8_UNIFORM_01 returns a unit pseudorandom R8.
  //
  //  Discussion:
  //
  //    This routine implements the recursion
  //
  //      seed = 16807 * seed mod ( 2^31 - 1 )
  //      unif = seed / ( 2^31 - 1 )
  //
  //    The integer arithmetic never requires more than 32 bits,
  //    including a sign bit.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    11 August 2004
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Paul Bratley, Bennett Fox, Linus Schrage,
  //    A Guide to Simulation,
  //    Springer Verlag, pages 201-202, 1983.
  //
  //    Bennett Fox,
  //    Algorithm 647:
  //    Implementation and Relative Efficiency of Quasirandom
  //    Sequence Generators,
  //    ACM Transactions on Mathematical Software,
  //    Volume 12, Number 4, pages 362-376, 1986.
  //
  //  Parameters:
  //
  //    Input/output, int &SEED, the "seed" value.  Normally, this
  //    value should not be 0.  On output, SEED has been updated.
  //
  //    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
  //    0 and 1.
  //
{
  int k;
  double r;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8poly_value_horner ( int m, double c[], double x )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
  //
  //  Discussion:
  //
  //    The polynomial
  //
  //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
  //
  //    is to be evaluated at the value X.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    02 January 2015
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int M, the degree of the polynomial.
  //
  //    Input, double C[M+1], the coefficients of the polynomial.
  //    A[0] is the constant term.
  //
  //    Input, double X, the point at which the polynomial is to be evaluated.
  //
  //    Output, double R8POLY_VALUE_HORNER, the value of the polynomial at X.
  //
{
  int i;
  double value;

  value = c[m];

  for ( i = m - 1; 0 <= i; i-- )
  {
    value = value * x + c[i];
  }

  return value;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
  //
  //  Discussion:
  //
  //    An R8VEC is a vector of R8's.
  //
  //    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
  //
  //    In other words, the interval is divided into N-1 even subintervals,
  //    and the endpoints of intervals are used as the points.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    29 March 2011
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the number of entries in the vector.
  //
  //    Input, double A_FIRST, A_LAST, the first and last entries.
  //
  //    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
  //
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first
                 + ( double ) (         i ) * a_last )
      / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double r8vec_max ( int n, double *dvec )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8VEC_MAX returns the value of the maximum element in an R8VEC.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    15 October 1998
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the number of entries in the array.
  //
  //    Input, double *RVEC, a pointer to the first entry of the array.
  //
  //    Output, double R8VEC_MAX, the value of the maximum element.  This
  //    is set to 0.0 if N <= 0.
  //
{
  int i;
  double rmax;
  double *r8vec_pointer;

  if ( n <= 0 )
  {
    return 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      rmax = *dvec;
      r8vec_pointer = dvec;
    }
    else
    {
      r8vec_pointer++;
      if ( rmax < *r8vec_pointer )
      {
        rmax = *r8vec_pointer;
      }
    }
  }

  return rmax;

}
//****************************************************************************80

double r8vec_mean ( int n, double x[] )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8VEC_MEAN returns the mean of an R8VEC.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    01 May 1999
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the number of entries in the vector.
  //
  //    Input, double X[N], the vector whose mean is desired.
  //
  //    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
  //
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
//****************************************************************************80

double r8vec_min ( int n, double *dvec )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8VEC_MIN returns the value of the minimum element in an R8VEC.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    15 October 1998
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the number of entries in the array.
  //
  //    Input, double *RVEC, a pointer to the first entry of the array.
  //
  //    Output, double R8VEC_MIN, the value of the minimum element.  This
  //    is set to 0.0 if N <= 0.
  //
{
  int i;
  double rmin;
  double *r8vec_pointer;

  if ( n <= 0 )
  {
    return 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      rmin = *dvec;
      r8vec_pointer = dvec;
    }
    else
    {
      r8vec_pointer++;
      if ( *r8vec_pointer < rmin )
      {
        rmin = *r8vec_pointer;
      }
    }
  }

  return rmin;
}
//****************************************************************************80

double r8vec_variance ( int n, double x[] )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    R8VEC_VARIANCE returns the variance of an R8VEC.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    01 May 1999
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int N, the number of entries in the vector.
  //
  //    Input, double X[N], the vector whose variance is desired.
  //
  //    Output, double R8VEC_VARIANCE, the variance of the vector entries.
  //
{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
//****************************************************************************80
double truncated_normal_ab_cdf ( double x, double mu, double sigma, double a,
                                 double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_CDF evaluates the truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the CDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Output, double TRUNCATED_NORMAL_AB_CDF, the value of the CDF.
  //
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );
  xi_cdf = normal_01_cdf ( xi );

  cdf = ( xi_cdf - alpha_cdf ) / ( beta_cdf - alpha_cdf );

  return cdf;
}
//****************************************************************************80

void truncated_normal_ab_cdf_values ( int &n_data, double &mu, double &sigma,
                                      double &a, double &b, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_CDF_VALUES: values of the Truncated Normal CDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval [A,B].
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    13 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &A, &B, the lower and upper truncation limits.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0 };

  static double b_vec[N_MAX] = {
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0 };

  static double fx_vec[N_MAX] = {
    0.3371694242213513,
    0.3685009225506048,
    0.4006444233448185,
    0.4334107066903040,
    0.4665988676496338,
    0.5000000000000000,
    0.5334011323503662,
    0.5665892933096960,
    0.5993555766551815,
    0.6314990774493952,
    0.6628305757786487 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_ab_cdf_inv ( double cdf, double mu, double sigma, double a,
                                     double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_CDF_INV inverts the truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double CDF, the value of the CDF.
  //    0.0 <= CDF <= 1.0.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Output, double TRUNCATED_NORMAL_AB_CDF_INV, the corresponding argument.
  //
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double x;
  double xi;
  double xi_cdf;


  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  xi_cdf = ( beta_cdf - alpha_cdf ) * cdf + alpha_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_ab_mean ( double mu, double sigma, double a, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_MEAN returns the mean of the truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Output, double TRUNCATED_NORMAL_AB_MEAN, the mean of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double beta;
  double beta_cdf;
  double beta_pdf;
  double mean;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta  );

  alpha_pdf = normal_01_pdf ( alpha );
  beta_pdf = normal_01_pdf ( beta );

  mean = mu + sigma * ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf );

  return mean;
}
//****************************************************************************80

double truncated_normal_ab_moment ( int order, double mu, double sigma,
                                    double a, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    11 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Phoebus Dhrymes,
  //    Moments of Truncated Normal Distributions,
  //    May 2005.
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //    0.0 < S.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //    A < B.
  //
  //    Output, double TRUNCATED_NORMAL_AB_MOMENT, the moment of the PDF.
  //
{
  double a_cdf;
  double a_h;
  double a_pdf;
  double b_cdf;
  double b_h;
  double b_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;

  a_h = ( a - mu ) / sigma;
  a_pdf = normal_01_pdf ( a_h );
  a_cdf = normal_01_cdf ( a_h );


  b_h = ( b - mu ) / sigma;
  b_pdf = normal_01_pdf ( b_h );
  b_cdf = normal_01_cdf ( b_h );


  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf );
    }
    else
    {
      ir = ( double ) ( r - 1 ) * irm2
      - ( pow ( b_h, r - 1 ) * b_pdf - pow ( a_h, r - 1 ) * a_pdf )
      / ( b_cdf - a_cdf );
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r )
      * pow ( sigma, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
//****************************************************************************80

double truncated_normal_ab_pdf ( double x, double mu, double sigma, double a,
                                 double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_PDF evaluates the truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the PDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Output, double TRUNCATED_NORMAL_AB_PDF, the value of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );
  xi_pdf = normal_01_pdf ( xi );

  pdf = xi_pdf / ( beta_cdf - alpha_cdf ) / sigma;

  return pdf;
}
//****************************************************************************80

void truncated_normal_ab_pdf_values ( int &n_data, double &mu, double &sigma,
                                      double &a, double &b, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_PDF_VALUES: values of the Truncated Normal PDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval [A,B].
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    13 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &A, &B, the lower and upper truncation limits.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0 };

  static double b_vec[N_MAX] = {
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0 };

  static double fx_vec[N_MAX] = {
    0.01543301171801836,
    0.01588394472270638,
    0.01624375997031919,
    0.01650575046469259,
    0.01666496869385951,
    0.01671838200940538,
    0.01666496869385951,
    0.01650575046469259,
    0.01624375997031919,
    0.01588394472270638,
    0.01543301171801836 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_ab_sample ( double mu, double sigma, double a, double b,
                                    int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Input/output, int &SEED, a seed for the random number
  //    generator.
  //
  //    Output, double TRUNCATED_NORMAL_AB_SAMPLE, a sample of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  u = r8_uniform_01 ( seed );
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_ab_variance ( double mu, double sigma, double a,
                                      double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_AB_VARIANCE returns the variance of the truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, B, the lower and upper truncation limits.
  //
  //    Output, double TRUNCATED_NORMAL_AB_VARIANCE, the variance of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double beta;
  double beta_cdf;
  double beta_pdf;
  double variance;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_pdf = normal_01_pdf ( alpha );
  beta_pdf = normal_01_pdf ( beta );

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  variance = sigma * sigma * ( 1.0
                                 + ( alpha * alpha_pdf - beta * beta_pdf ) / ( beta_cdf - alpha_cdf )
                                 - pow ( ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf ), 2 ) );

                                 return variance;
}
//****************************************************************************80

double truncated_normal_a_cdf ( double x, double mu, double sigma, double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_CDF evaluates the lower truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the CDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_CDF, the value of the CDF.
  //
{
  double alpha;
  double alpha_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  xi_cdf = normal_01_cdf ( xi );

  cdf = ( xi_cdf - alpha_cdf ) / ( 1.0 - alpha_cdf );

  return cdf;
}
//****************************************************************************80

void truncated_normal_a_cdf_values ( int &n_data, double &mu, double &sigma,
                                     double &a, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_CDF_VALUES: values of the Lower Truncated Normal CDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval [A,+oo).
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &A, the lower truncation limit.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0 };

  static double fx_vec[N_MAX] = {
    0.3293202045481688,
    0.3599223134505957,
    0.3913175216041539,
    0.4233210140873113,
    0.4557365629792204,
    0.4883601253415709,
    0.5209836877039214,
    0.5533992365958304,
    0.5854027290789878,
    0.6167979372325460,
    0.6474000461349729 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_a_cdf_inv ( double cdf, double mu, double sigma,
                                    double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_CDF_INV inverts the lower truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double CDF, the value of the CDF.
  //    0.0 <= CDF <= 1.0.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_CDF_INV, the corresponding argument.
  //
{
  double alpha;
  double alpha_cdf;
  double x;
  double xi;
  double xi_cdf;



  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  xi_cdf = ( 1.0 - alpha_cdf ) * cdf + alpha_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_a_mean ( double mu, double sigma, double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_MEAN returns the mean of the lower truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_MEAN, the mean of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double mean;

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  alpha_pdf = normal_01_pdf ( alpha );

  mean = mu + sigma * alpha_pdf / ( 1.0 - alpha_cdf );

  return mean;
}
//****************************************************************************80

double truncated_normal_a_moment ( int order, double mu, double sigma,
                                   double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    11 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Phoebus Dhrymes,
  //    Moments of Truncated Normal Distributions,
  //    May 2005.
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_MOMENT, the moment of the PDF.
  //
{
  double moment;

  moment = r8_mop ( order )
    * truncated_normal_b_moment ( order, - mu, sigma, - a );

  return moment;
}
//****************************************************************************80

double truncated_normal_a_pdf ( double x, double mu, double sigma, double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_PDF evaluates the lower truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the PDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_PDF, the value of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  alpha = ( a - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  xi_pdf = normal_01_pdf ( xi );

  pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma;

  return pdf;
}
//****************************************************************************80

void truncated_normal_a_pdf_values ( int &n_data, double &mu, double &sigma,
                                     double &a, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_PDF_VALUES: values of the Lower Truncated Normal PDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval [A,+oo).
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &A, the lower truncation limit.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0,
    50.0 };

  static double fx_vec[N_MAX] = {
    0.01507373507401876,
    0.01551417047139894,
    0.01586560931024694,
    0.01612150073158793,
    0.01627701240029317,
    0.01632918226724295,
    0.01627701240029317,
    0.01612150073158793,
    0.01586560931024694,
    0.01551417047139894,
    0.01507373507401876 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_a_sample ( double mu, double sigma, double a,
                                   int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_SAMPLE samples the lower truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Input/output, int &SEED, a seed for the random number
  //    generator.
  //
  //    Output, double TRUNCATED_NORMAL_A_SAMPLE, a sample of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  u = r8_uniform_01 ( seed );
  xi_cdf = alpha_cdf + u * ( 1.0 - alpha_cdf );
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_a_variance ( double mu, double sigma, double a )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_A_VARIANCE: variance of the lower truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double A, the lower truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_A_VARIANCE, the variance of the PDF.
  //
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double variance;

  alpha = ( a - mu ) / sigma;

  alpha_pdf = normal_01_pdf ( alpha );

  alpha_cdf = normal_01_cdf ( alpha );

  variance = sigma * sigma * ( 1.0
                                 + alpha * alpha_pdf / ( 1.0 - alpha_cdf )
                                 - pow ( alpha_pdf / ( 1.0 - alpha_cdf ), 2 ) );

                                 return variance;
}
//****************************************************************************80

double truncated_normal_b_cdf ( double x, double mu, double sigma, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_CDF evaluates the upper truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the CDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_CDF, the value of the CDF.
  //
{
  double beta;
  double beta_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  beta = ( b - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );
  xi_cdf = normal_01_cdf ( xi );

  cdf = xi_cdf / beta_cdf;

  return cdf;
}
//****************************************************************************80

void truncated_normal_b_cdf_values ( int &n_data, double &mu, double &sigma,
                                     double &b, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_CDF_VALUES: values of the upper Truncated Normal CDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval (-oo,B].
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &B, the upper truncation limit.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double b_vec[N_MAX] = {
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0 };

  static double fx_vec[N_MAX] = {
    0.3525999538650271,
    0.3832020627674540,
    0.4145972709210122,
    0.4466007634041696,
    0.4790163122960786,
    0.5116398746584291,
    0.5442634370207796,
    0.5766789859126887,
    0.6086824783958461,
    0.6400776865494043,
    0.6706797954518312 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_b_cdf_inv ( double cdf, double mu, double sigma,
                                    double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_CDF_INV inverts the upper truncated Normal CDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double CDF, the value of the CDF.
  //    0.0 <= CDF <= 1.0.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_CDF_INV, the corresponding argument.
  //
{
  double beta;
  double beta_cdf;
  double x;
  double xi;
  double xi_cdf;



  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  xi_cdf = beta_cdf * cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_b_mean ( double mu, double sigma, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_MEAN returns the mean of the upper truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_MEAN, the mean of the PDF.
  //
{
  double beta;
  double beta_cdf;
  double beta_pdf;
  double mean;

  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  beta_pdf = normal_01_pdf ( beta );

  mean = mu - sigma * beta_pdf / beta_cdf;

  return mean;
}
//****************************************************************************80

double truncated_normal_b_moment ( int order, double mu, double sigma, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    11 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Phoebus Dhrymes,
  //    Moments of Truncated Normal Distributions,
  //    May 2005.
  //
  //  Parameters:
  //
  //    Input, int ORDER, the order of the moment.
  //    0 <= ORDER.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_MOMENT, the moment of the PDF.
  //
{
  double f;
  double h;
  double h_cdf;
  double h_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;



  h = ( b - mu ) / sigma;
  h_pdf = normal_01_pdf ( h );
  h_cdf = normal_01_cdf ( h );



  f = h_pdf / h_cdf;

  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - f;
    }
    else
    {
      ir = - pow ( h, r - 1 ) * f + ( double ) ( r - 1 ) * irm2;
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r )
      * pow ( sigma, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
//****************************************************************************80

double truncated_normal_b_pdf ( double x, double mu, double sigma, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_PDF evaluates the upper truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double X, the argument of the PDF.
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_PDF, the value of the PDF.
  //
{
  double beta;
  double beta_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  beta = ( b - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );
  xi_pdf = normal_01_pdf ( xi );

  pdf = xi_pdf / beta_cdf / sigma;

  return pdf;
}
//****************************************************************************80

void truncated_normal_b_pdf_values ( int &n_data, double &mu, double &sigma,
                                     double &b, double &x, double &fx )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_PDF_VALUES: values of the Upper Truncated Normal PDF.
  //
  //  Discussion:
  //
  //    The Normal distribution, with mean Mu and standard deviation Sigma,
  //    is truncated to the interval (-oo,B].
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 September 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
  //    first call.  On each call, the routine increments N_DATA by 1, and
  //    returns the corresponding data; when there is no more data, the
  //    output value of N_DATA will be 0 again.
  //
  //    Output, double &MU, the mean of the distribution.
  //
  //    Output, double &SIGMA, the standard deviation of the distribution.
  //
  //    Output, double &B, the upper truncation limit.
  //
  //    Output, double &X, the argument of the function.
  //
  //    Output, double &FX, the value of the function.
  //
{
# define N_MAX 11

  static double b_vec[N_MAX] = {
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0,
    150.0 };

  static double fx_vec[N_MAX] = {
    0.01507373507401876,
    0.01551417047139894,
    0.01586560931024694,
    0.01612150073158793,
    0.01627701240029317,
    0.01632918226724295,
    0.01627701240029317,
    0.01612150073158793,
    0.01586560931024694,
    0.01551417047139894,
    0.01507373507401876 };

  static double mu_vec[N_MAX] = {
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0,
    100.0 };

  static double sigma_vec[N_MAX] = {
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0,
    25.0  };

  static double x_vec[N_MAX] = {
    90.0,
    92.0,
    94.0,
    96.0,
    98.0,
    100.0,
    102.0,
    104.0,
    106.0,
    108.0,
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_b_sample ( double mu, double sigma, double b, int &seed )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_SAMPLE samples the upper truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    21 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Input/output, int &SEED, a seed for the random number
  //    generator.
  //
  //    Output, double TRUNCATED_NORMAL_B_SAMPLE, a sample of the PDF.
  //
{
  double beta;
  double beta_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  u = r8_uniform_01 ( seed );
  xi_cdf = u * beta_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_b_variance ( double mu, double sigma, double b )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TRUNCATED_NORMAL_B_VARIANCE: variance of the upper truncated Normal PDF.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 August 2013
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, double MU, SIGMA, the mean and standard deviation of the
  //    parent Normal distribution.
  //
  //    Input, double B, the upper truncation limit.
  //
  //    Output, double TRUNCATED_NORMAL_B_VARIANCE, the variance of the PDF.
  //
{
  double beta;
  double beta_cdf;
  double beta_pdf;
  double variance;

  beta = ( b - mu ) / sigma;

  beta_pdf = normal_01_pdf ( beta );

  beta_cdf = normal_01_cdf ( beta );

  variance = sigma * sigma * ( 1.0
                                 - beta * beta_pdf / beta_cdf
                                 - pow ( beta_pdf / beta_cdf, 2 ) );

                                 return variance;
}












//****************************************************************************80
arma::rowvec nonnan_ind(arma::rowvec x){
  int n=x.n_cols;
  arma::rowvec ans(n);
  int flag = 0;
  for(int i=0;i<n; i++){
    if(x(i)==x(i)){
      ans(flag) = i;
      flag ++;
    }
  }
  ans.reshape(1,flag);
  return (ans);
}

arma::mat rowsome(arma::mat x,
                  arma::rowvec ind){
  int C = x.n_cols;
  int R = ind.n_cols;
  arma::mat ans(R,C);
  for(int i=0; i<R; i++){
    ans.row(i) = x.row(ind(i));
  }
  return(ans);
}

arma::rowvec ton(int n){
  //ton(4) should give (0,1,2,3)
  arma::rowvec ans(n);
  for(int i=0; i<n; i++){
    ans(i) = i;
  }
  return ans;
}

arma::mat colsome(arma::mat x,
                  arma::rowvec ind){
  int R = x.n_rows;
  int C = ind.n_cols;
  arma::mat ans(R,C);
  for(int i=0; i<C; i++){
    ans.col(i) = x.col(ind(i));
  }
  return(ans);
}

arma::mat subsome(arma::mat x,
                  arma::rowvec ind){
  int n = ind.n_cols;
  int m = x.n_rows;
  arma::mat temp(m, n);
  arma::mat ans(n,n);
  temp = rowsome(x, ind);
  ans = colsome(temp,ind);
  return (ans);
}

double vectornorm(arma::rowvec x){
  double ans=0;
  ans = sum(x%x);
  return ans;
}

double dmvnrm_arma0(arma::rowvec x,
                    arma::rowvec mean,
                    double sigma_square,
                    bool logd = false) {

  double out = vectornorm(x-mean);
  out = out/(-2 * sigma_square);
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

double dmvnrm_arma(arma::rowvec x,
                   arma::rowvec mean,
                   double sigma_square,
                   bool logd = false) {
  arma::rowvec non_ind = nonnan_ind(x);
  arma::rowvec xnew = colsome(x,non_ind);
  arma::rowvec meannew = colsome(mean, non_ind);

  if(xnew.n_cols == 0){
    if(logd){
      return (0);
    }
    else{
      return (1);
    }
  }

  return(dmvnrm_arma0(xnew,meannew,sigma_square,logd));
}

int rmunoim(arma::rowvec probs) {
  //this function returns a int, not a vector
  //rmultinom(1,1,probs)
  int k = probs.n_cols;
  IntegerVector ans(k);
  int r=0;
  rmultinom(1, probs.begin(), k, ans.begin());
  for(int i=0; i<k; i++){
    if(ans(i)==1){
      r=i+1;
    }
  }
  return r;
}

int sampleint(int n){
  arma::rowvec beta(n);
  arma::mat temp(1,1);
  temp = 1.0/n;
  beta = arma::repmat(temp,1,n);
  int ans = rmunoim(beta);
  return ans;
}

double kernelC(arma::rowvec x, arma::rowvec y, double theta, double tau){
  double ans;
  ans = exp(-vectornorm(x-y)/(theta*theta));
  return ans;
}

arma::mat eye(int n){
  arma:: mat ans(n,n);
  ans.eye();
  return (ans);
}

arma::mat updateC(arma::mat Z,double theta, double tau){
  int K=Z.n_rows;
  int L = Z.n_cols;
  int j;
  arma::rowvec q(K);
  arma::mat C = eye(K);
  arma::rowvec mean(L);
  arma::mat Temp(K,K);
  mean.fill(0);
  C = C*0.00001;

  for(int i=0;i<K;i++){
    q(i) = dmvnrm_arma0(Z.row(i),mean,tau*tau);
    j = K-1;
    while(j>=i){
      C(i,j) += kernelC(Z.row(i),Z.row(j),theta,tau);
      j--;
    }
    while(j>=0){
      C(i,j) = C(j,i);
      j--;
    }
  }
  Temp = q.t()*q;
  return (C % Temp);
  //return C;
}


arma::rowvec getind(arma::rowvec x, double c){
  int n = x.n_cols;
  arma::rowvec ans(n);
  int flag = 0;
  for(int i=0; i<n; i++){
    if(x(i)==c){
      ans(flag) = i;
      flag++;
    }
  }
  ans.reshape(1,flag);
  return(ans);
}

arma::mat matrix(arma::rowvec x, int n){
  int m = x.n_cols;
  arma::mat ans(n,m);
  for(int i=0; i<n; i++){
    ans.row(i)=x;
  }
  return(ans);
}

arma::rowvec getgamma(arma::cube x, int i, int j){
  int n = x.n_cols;
  int m = x.n_rows;
  arma::mat temp(m,n);
  arma::rowvec ans(n);
  temp = x.slice(i);
  ans = temp.row(j);
  return(ans);
}

arma::mat minusvalue(arma::mat x, int j, arma::rowvec ind, arma::rowvec y){
  //used for Y.theta[i,indj]<-Y[i,indj]-xbeta
  int n=ind.n_cols;
  for(int i=0; i<n;i++){
    x(j,ind(i)) -= y(i);
  }
  return x;
}

double findT(arma::rowvec x){
  double ans= -1;
  int n = x.n_cols;
  for(int i=0; i<n; i++){
    if(x(i)>0){
      ans = i;
    }
  }
  return ans;
}

arma::mat removeii(arma::mat A, int a, int b){
  //a,b start from 1
  //remove ath row and bth column
  int row=A.n_rows;
  int col=A.n_cols;
  arma::mat ans(row-1, col-1);
  arma::mat temp (row-1, col);
  for(int i=0;i<row-1;i++){
    if(i<a-1){
      temp.row(i) = A.row(i);
    }
    else{
      temp.row(i) = A.row(i+1);
    }
  }
  for(int j=0; j<col-1; j++){
    if(j<b-1){
      ans.col(j) = temp.col(j);
    }
    else{
      ans.col(j) = temp.col(j+1);
    }
  }


  return ans;
}

arma::rowvec removei(arma::rowvec A, int a){
  //need use i+1, since a starts from 1
  int k = A.n_cols;
  arma::rowvec ans(k-1);
  for (int i=0; i<k-1;i++){
    if(i<a-1){
      ans(i) = A(i);
    }
    else{
      ans(i) = A(i+1);
    }
  }
  return ans;
}

arma::vec cDmu(double x, arma::vec y){
  arma::rowvec ans(y.n_rows);
  ans = x * (y.t());
  return(ans.t());
}

arma::rowvec ind(arma::rowvec x, arma::rowvec index){
  int n = accu(index);
  arma::rowvec ans(n);
  int flag = 0;
  int i = 0;
  while(flag<n){
    if(index(i)){
      ans(flag) = x(i);
      flag ++;
    }
    i++;
  }
  return(ans);
}

arma::mat choleskyDecompose(arma::mat A){
  double sum;
  int size=A.n_cols ;
  arma::mat tmp(size, size);
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      if(i<j){
        tmp(i,j) = 0.0;
        continue;
      }
      sum = 0.0;
      for(int k=0;k<=j-1;k++){
        sum = sum+tmp(i,k) *tmp(j,k);
      }
      if(j<i){
        tmp(i,j) = (A(i,j) - sum) / tmp(j,j);
        //temp[i][j]=(elements[i][j]-sum)/temp[j][j];
      }
      else if(j==i){
        tmp(i,j) = sqrt(A(i,j)-sum);
        //temp[i][j]=sqrt(elements[i][j]-sum);
      }
    }
  }
  return tmp;
}

arma::mat inverse_lower(arma::mat A){
  int size=A.n_cols;
  arma::mat tmp = eye(size);

  for(int i=0; i<size; i++){
    tmp.row(i) = tmp.row(i)/A(i,i);
    A.row(i) = A.row(i)/A(i,i);

    for(int j=i+1; j<size; j++){
      tmp.row(j) = tmp.row(j) - tmp.row(i) * A(j,i);
      A.row(j) = A.row(j) - A.row(i);
    }

  }
  return tmp;
}

arma::mat inverse_via_chol(arma::mat A){
  int size=A.n_cols;
  arma::mat chol = choleskyDecompose(A);
  arma::mat tmp = inverse_lower(chol);
  arma::mat ans = tmp.t() * tmp;
  return ans;
}


List updateBeta(arma::mat X,
                arma::mat Y,
                arma::mat Z,
                arma::mat delta,
                arma::mat Beta,
                arma::cube Gamma,
                arma::rowvec E,
                arma::mat R,
                double S,
                arma::rowvec Ds,
                arma::mat mustar,
                arma::mat mu,
                double sigma,
                arma::rowvec c,
                arma::mat C,
                arma::mat step,
                arma::mat runif,
                int n, int m, int T0, int p, int q, arma::mat D,
                double theta, double tau){
  List ans;
  arma::rowvec Betanew(p);
  arma::vec mustarnewj(T0);
  arma::vec munewj(m);
  arma::mat mustarnew(n,T0);
  arma::mat munew(n,m);
  arma::mat Ybeta(n,m);
  arma::rowvec indj(m);
  arma::mat Zj(m,q);
  arma::rowvec gammaij(q);
  arma::rowvec betai(p);
  arma::rowvec indi(n);
  arma::vec beta0(m);
  arma::mat Cminus(m-1,m-1);
  arma::rowvec b(S-1);
  arma::mat temp(1,1);
  arma::rowvec bnew(S-1);
  int l;
  double likelinow, likelinew, fnow, fnew,va;
  munew = mu;
  mustarnew = mustar;
  Ybeta = Y;
  for(int i=0; i<n; i++){
    for(int j=1; j<=Ds(E(i)-1);j++){
      l=getind(R.row(E(i)-1),j).n_cols;
      indj.set_size(1,l);
      indj = getind(R.row(E(i)-1),j);
      Zj.set_size(l,q);
      Zj = rowsome(Z,indj);
      gammaij = Gamma.slice(E(i)-1).row(j-1);
      Ybeta = minusvalue(Ybeta,i,indj, (Zj*gammaij.t()).t());
    }
  }
  //now update beta one by one row
  for(int i=0; i<S;i++){
    likelinow = 0;
    betai = Beta.row(i);
    l=getind(E,i+1).n_cols;
    indi.set_size(1,l);
    indi = getind(E,i+1);
    for(int j=0; j<l; j++){
      //j in R is indi(j) now
      beta0 = matrix(X.row(indi(j)),m) * betai.t();
      likelinow += dmvnrm_arma(Ybeta.row(indi(j)),  beta0.t(), sigma , true);
      for(int t=0; t<T0; t++){
        va = min(5.0, mustar(indi(j),t));
        va = max(va, -5.0);
        likelinow += R::pnorm(va,0,1,true, true)* delta(indi(j),t) +
          log(1-R::pnorm(va,0,1,true, false))*(1-delta(indi(j),t));
      }
    }
    Cminus = removeii(C,i+1,i+1);
    b = removei(C.row(i),i+1);
    temp = b*inverse_via_chol(Cminus)*b.t();
    fnow = kernelC(Beta.row(i),Beta.row(i),theta,tau) -  temp(0,0);
    //update Beta one by one
    for(int j=0; j<p;j++){
      Betanew = Beta.row(i);
      Betanew(j) = Betanew(j) + step(i,j);
      likelinew = 0;
      for(int k=0; k<l;k++){
        //k in R is now indi(k)
        beta0 = matrix(X.row(indi(k)),m)*Betanew.t();

        likelinew += dmvnrm_arma(Ybeta.row(indi(k)),  beta0.t(),  sigma , true);
        munewj = mu.row(indi(k)).t() + matrix(X.row(indi(k)),m) * (Betanew.t()-Beta.row(i).t());
        mustarnewj = mustar.row(indi(k)).t() + cDmu(c(1),D*(munewj)) - cDmu(c(1),D*mu.row(indi(k)).t());
        mustarnew.row(indi(k)) = mustarnewj.t();
        munew.row(indi(k)) = munewj.t();
        for(int t=0; t<T0; t++){
          va = min(5.0, mustarnewj(t));
          va = max(va, -5.0);
          likelinew += R::pnorm(va,0,1,true,true)*delta(indi(k),t) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(indi(k),t));
        }
      }
      //compute prior
      for(int s=0; s<S; s++){
        if(s<i){
          bnew(s) = kernelC(Beta.row(s), Betanew, theta, tau);
        }
        else if(s>i){
          bnew(s-1) = kernelC(Beta.row(s), Betanew,theta,tau);
        }
      }
      temp = bnew*inverse_via_chol(Cminus)*bnew.t();
      fnew = kernelC(Betanew, Betanew,theta,tau) - temp(0,0);
      if(l==0){
        likelinew = 0;
        likelinow = 0;
      }
      //cout<<"i="<<i+1<<"j="<<j+1<<"likelinow="<<likelinow<<", likelinew="<<likelinew<<endl;
      if((likelinew - likelinow + log(fnew/fnow))>log(runif(i,j))){
        b = bnew;
        fnow = fnew;
        //update C
        for(int s=0; s<S; s++){
          if(s<i){
            C(i,s) = bnew(s);
            C(s,i) = bnew(s);
          }
          else if(s>i){
            C(i,s) = bnew(s-1);
            C(s,i) = bnew(s-1);
          }
          else{
            C(i,i) = kernelC(Betanew,Betanew,theta, tau);
          }
        }
        likelinow = likelinew;
        Beta.row(i) = Betanew;
        mustar = mustarnew;
        mu = munew;
      }
    }
  }
  ans["C"] = C;
  ans["Beta"] = Beta;
  ans["mu"]= mu;
  ans["mustar"] = mustar;
  ans["Ybeta"] = Ybeta;
  return (ans);

}


List updateBetanomiss(arma::mat X,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::mat Beta,
                      arma::cube Gamma,
                      arma::rowvec E,
                      arma::mat R,
                      double S,
                      arma::rowvec Ds,
                      arma::mat mu,
                      double sigma,
                      arma::mat C,
                      arma::mat step,
                      arma::mat runif,
                      int n, int m, int T0, int p, int q, arma::mat D,
                      double theta, double tau){
  List ans;
  arma::rowvec Betanew(p);
  arma::vec munewj(m);
  arma::mat munew(n,m);
  arma::mat Ybeta(n,m);
  arma::rowvec indj(m);
  arma::mat Zj(m,q);
  arma::rowvec gammaij(q);
  arma::rowvec betai(p);
  arma::rowvec indi(n);
  arma::vec beta0(m);
  arma::mat Cminus(m-1,m-1);
  arma::rowvec b(S-1);
  arma::mat temp(1,1);
  arma::rowvec bnew(S-1);
  int l;
  double likelinow, likelinew, fnow, fnew;
  munew = mu;
  Ybeta = Y;
  for(int i=0; i<n; i++){
    for(int j=1; j<=Ds(E(i)-1);j++){
      l=getind(R.row(E(i)-1),j).n_cols;
      indj.set_size(1,l);
      indj = getind(R.row(E(i)-1),j);
      Zj.set_size(l,q);
      Zj = rowsome(Z,indj);
      gammaij = Gamma.slice(E(i)-1).row(j-1);
      Ybeta = minusvalue(Ybeta,i,indj, (Zj*gammaij.t()).t());
    }
  }
  //now update beta one by one row
  for(int i=0; i<S;i++){
    likelinow = 0;
    betai = Beta.row(i);
    l=getind(E,i+1).n_cols;
    indi.set_size(1,l);
    indi = getind(E,i+1);
    for(int j=0; j<l; j++){
      //j in R is indi(j) now
      beta0 = matrix(X.row(indi(j)),m) * betai.t();
      likelinow += dmvnrm_arma(Ybeta.row(indi(j)),  beta0.t(), sigma , true);
    }
    Cminus = removeii(C,i+1,i+1);
    b = removei(C.row(i),i+1);
    temp = b*inverse_via_chol(Cminus)*b.t();
    fnow = kernelC(Beta.row(i),Beta.row(i),theta,tau) -  temp(0,0);
    //update Beta one by one
    for(int j=0; j<p;j++){
      Betanew = Beta.row(i);
      Betanew(j) = Betanew(j) + step(i,j);
      likelinew = 0;
      for(int k=0; k<l;k++){
        //k in R is now indi(k)
        beta0 = matrix(X.row(indi(k)),m)*Betanew.t();

        likelinew += dmvnrm_arma(Ybeta.row(indi(k)),  beta0.t(),  sigma , true);
        munewj = mu.row(indi(k)).t() + matrix(X.row(indi(k)),m) * (Betanew.t()-Beta.row(i).t());
        munew.row(indi(k)) = munewj.t();
      }
      //compute prior
      for(int s=0; s<S; s++){
        if(s<i){
          bnew(s) = kernelC(Beta.row(s), Betanew, theta, tau);
        }
        else if(s>i){
          bnew(s-1) = kernelC(Beta.row(s), Betanew,theta,tau);
        }
      }
      temp = bnew*inverse_via_chol(Cminus)*bnew.t();
      fnew = kernelC(Betanew, Betanew,theta,tau) - temp(0,0);
      if(l==0){
        likelinew = 0;
        likelinow = 0;
      }
      //cout<<"i="<<i+1<<"j="<<j+1<<"likelinow="<<likelinow<<", likelinew="<<likelinew<<endl;
      if((likelinew - likelinow + log(fnew/fnow))>log(runif(i,j))){
        b = bnew;
        fnow = fnew;
        //update C
        for(int s=0; s<S; s++){
          if(s<i){
            C(i,s) = bnew(s);
            C(s,i) = bnew(s);
          }
          else if(s>i){
            C(i,s) = bnew(s-1);
            C(s,i) = bnew(s-1);
          }
          else{
            C(i,i) = kernelC(Betanew,Betanew,theta, tau);
          }
        }
        likelinow = likelinew;
        Beta.row(i) = Betanew;
        mu = munew;
      }
    }
  }
  ans["C"] = C;
  ans["Beta"] = Beta;
  ans["mu"]= mu;
  return (ans);

}

int Nuni(arma::rowvec x){
  arma::rowvec a = unique(x);
  int n = a.n_cols;
  return (n);
}

arma::mat getGammai(arma::cube Gamma, int j, int ds, int q){
  arma::mat ans(ds, q);
  for(int i=0; i<ds;i++){
    ans.row(i) = getgamma(Gamma, j,i);
  }
  return(ans);
}

arma::mat putvalue(arma::mat x, int k, arma::rowvec indi, arma::rowvec value){
  int n = indi.n_cols;
  for(int i=0; i<n; i++){
    x(k,indi(i)) = value(i);
  }
  return (x);
}

arma::rowvec computemu(arma::mat mu, int k, arma::rowvec indi, arma::mat Zi, arma::rowvec gamma, arma::rowvec gammanew){
  int n = indi.n_cols;
  arma::rowvec ans(n);
  arma::vec temp(n);
  ans = colsome(mu.row(k), indi);
  temp = Zi *(gammanew.t()-gamma.t());
  ans += temp.t();
  return (ans);
}

arma::rowvec computemustar(arma::mat mustar, int k, arma::rowvec indT, arma::rowvec c,
                           arma::mat D, arma::rowvec muk, arma::rowvec muknew, int T0){
  int n = indT.n_cols;
  arma::vec temp(T0);
  arma::rowvec ans(n);
  ans = colsome(mustar.row(k), indT);
  temp = cDmu(c(1),D*muknew.t()) - cDmu(c(1), D*muk.t());
  ans += colsome(temp.t(), indT);
  return(ans);
}

arma::cube putGamma(arma::cube Gamma, arma::mat Gammai, int i, int q){
  arma::cube ans = Gamma;
  int R = Gamma.n_rows;
  int n = Gammai.n_rows;
  arma::mat temp(R,q);
  temp = ans.slice(i);
  for(int j=0; j<n; j++){
    temp.row(j) = Gammai.row(j);
  }
  ans.slice(i) = temp;
  return(ans);
}

arma::rowvec ind1(arma::rowvec x, arma::rowvec indi, arma::rowvec index){
  //indi starts from 0
  int l = indi.n_cols;
  arma::rowvec ans(l);
  ans.fill(0);
  int flag = 0;
  for(int i=0; i<l; i++){
    if(index(indi(i))){
      ans(flag) = x(indi(i));
      flag++;
    }
  }
  ans.reshape(1,flag);
  return(ans);
}

arma::mat removeiii(arma::mat x, int a){
  //a starts from 1
  int R = x.n_rows-1;
  int C = x.n_cols;
  arma::mat ans(R,C);
  for(int i=0; i<R; i++){
    if(i<a-1){
      ans.row(i) = x.row(i);
    }
    else{
      ans.row(i) = x.row(i+1);
    }
  }
  return (ans);
}


List updateGamma(arma::mat X,
                 arma::mat Y,
                 arma::mat Z,
                 arma::mat delta,
                 arma::mat Beta,
                 arma::cube Gamma,
                 arma::rowvec E,
                 arma::mat R,
                 double S,
                 arma::rowvec Ds,
                 arma::mat mu,
                 arma::mat mustar,
                 double sigma,
                 arma::rowvec c,
                 arma::cube step,
                 arma::cube runif,
                 int n, int m, int T0, int p, int q, arma::mat D,
                 double theta, double tau){
  List ans;
  arma::mat mustarnew(n, T0);
  arma::mat munew(n, m);
  arma::rowvec inds(n);
  arma::rowvec indi(m);
  arma::rowvec indT(T0);
  arma::mat Gammai(max(Ds),q);
  mustarnew = mustar;
  arma::rowvec temp(T0);
  munew =mu;
  arma::mat C(10,10);
  arma::mat Cminus(9,9);
  arma::mat temp1(1,1);
  arma::rowvec b(9);
  arma::rowvec temp2(10);
  arma::rowvec bnew(9);
  arma::vec beta0(10);
  arma::rowvec Ygamma(m);
  arma::rowvec gammanew(q);
  //arma::cube L(max(Ds),q,S);
  //arma::cube P(max(Ds),q,S);
  int lens, ns, nt;
  double fnow, fnew, likelinow, likelinew, xbeta, va;
  for(int i=0; i<S; i++){
    //cout<<i<<endl;
    lens=getind(E,i+1).n_cols;
    inds.set_size(1,lens);
    inds = getind(E,i+1);
    Gammai.set_size(Ds(i),q);
    Gammai = getGammai(Gamma,i,Ds(i),q);
    C.set_size(Ds(i), Ds(i));
    C = updateC(Gammai,theta,tau);
    if(Ds(i)>1){
      Cminus.set_size(Ds(i)-1, Ds(i)-1);
      b.set_size(1,Ds(i)-1);
      temp2.set_size(1,Ds(i));
      bnew.set_size(1,Ds(i)-1);
    }
    for(int j=0; j<Ds(i); j++){
      //cout<<j<<endl;
      ns = getind(R.row(i),j+1).n_cols;
      indi.set_size(1,ns);
      indi = getind(R.row(i),j+1);
      beta0.set_size(ns,1);
      temp.set_size(1,ns);
      for(int k = 0; k<ns; k++){
        //k in R is indi(k) now
        temp(k) =  findT(D.col(indi(k)).t());
      }
      nt = Nuni(temp);
      indT.set_size(1,nt);
      indT = unique(temp);
      if(Ds(i)==1){
        //need to change is we change the kernel function
        fnow = 1;
      }
      else{
        Cminus = removeii(C,j+1,j+1);
        b = removei(C.row(j),j+1);
        temp1 = b*inverse_via_chol(Cminus)*b.t();
        fnow = kernelC(Gammai.row(j),Gammai.row(j),theta,tau) - temp1(0,0);
      }
      likelinow = 0;
      for(int k = 0; k<lens;k++){
        //k in R is now inds(k)
        xbeta = sum(X.row(inds(k))%Beta.row(i));
        Ygamma = Y.row(inds(k)) -  xbeta;
        beta0 = rowsome(Z,indi)* Gammai.row(j).t();
        likelinow += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);
        for(int t=0; t<nt; t++){
          //t in R is now indT(t)
          va = min(5.0, mustar(inds(k),indT(t)));
          va = max(va, -5.0);
          likelinow += R::pnorm(va,0,1,true, true)*delta(inds(k),indT(t)) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(inds(k),indT(t)));
        }
      }
      //update Gamma one element by one element
      for(int l=0; l<q; l++){
        gammanew = Gammai.row(j);
        gammanew(l) = gammanew(l) + step(j,l,i);
        likelinew = 0;
        for(int k=0; k<lens; k++){
          //k in R is now inds(k)
          xbeta = sum(X.row(inds(k))%Beta.row(i));
          Ygamma = Y.row(inds(k)) -  xbeta;
          beta0 = rowsome(Z,indi)* gammanew.t();
          likelinew += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);
          munew = putvalue(munew,inds(k),indi,
                           computemu(mu,inds(k),indi, rowsome(Z,indi), Gammai.row(j), gammanew));
          mustarnew = putvalue(mustarnew,inds(k),indT,
                               computemustar(mustar, inds(k),indT, c, D, mu.row(inds(k)), munew.row(inds(k)), T0));
          for(int t=0; t<nt; t++){
            //t in R is now indT(t)
            va = min(5.0, mustarnew(inds(k),indT(t)));
            va = max(va, -5.0);
            likelinew += R::pnorm(va,0,1,true, true)*delta(inds(k),indT(t)) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(inds(k),indT(t)));

          }

        }
        if(Ds(i)==1){
          fnew = 1;
          //need to change is we change the kernel function
        }
        else{
          for(int k=0; k<Ds(i);k++){
            temp2(k) = kernelC(Gammai.row(k), gammanew,theta,tau);
          }
          bnew = removei(temp2, j+1);
          temp1 = bnew*inverse_via_chol(Cminus)*bnew.t();
          fnew = kernelC(gammanew,gammanew,theta, tau) - temp1(0,0);
        }
        //L(j,l,i) = likelinew - likelinow;
        //P(j,l,i) = log(fnew/fnow);
        if((lens==0) || (ns==0)){
          likelinew = 0;
          likelinow = 0;
        }
        if((log(fnew/fnow) + likelinew - likelinow)>log(runif(j,l,i))){
          if(Ds(i)>1){
            fnow = fnew;
            C.row(j) = temp2;
            C.col(j) = temp2.t();
            C(j,j) = kernelC(gammanew,gammanew,theta,tau);
          }
          mu = munew;
          mustar = mustarnew;
          likelinow = likelinew;
          Gammai.row(j) = gammanew;
        }
        else{
          munew = mu;
          mustarnew = mustar;
        }
      }
    }
    //Gamma[,,i] = Gammai
    Gamma = putGamma(Gamma, Gammai, i, q);
  }
  ans["mu"]= mu;
  ans["mustar"] = mustar;
  ans["Gamma"]= Gamma;
  //ans["L"] = L;
  //ans["P"] = P;
  return(ans);
}

List updateGammanomiss(arma::mat X,
                       arma::mat Y,
                       arma::mat Z,
                       arma::mat delta,
                       arma::mat Beta,
                       arma::cube Gamma,
                       arma::rowvec E,
                       arma::mat R,
                       double S,
                       arma::rowvec Ds,
                       arma::mat mu,
                       double sigma,
                       arma::cube step,
                       arma::cube runif,
                       int n, int m, int T0, int p, int q, arma::mat D,
                       double theta, double tau){
  List ans;
  arma::mat mustarnew(n, T0);
  arma::mat munew(n, m);
  arma::rowvec inds(n);
  arma::rowvec indi(m);
  arma::rowvec indT(T0);
  arma::mat Gammai(max(Ds),q);
  arma::rowvec temp(T0);
  munew =mu;
  arma::mat C(10,10);
  arma::mat Cminus(9,9);
  arma::mat temp1(1,1);
  arma::rowvec b(9);
  arma::rowvec temp2(10);
  arma::rowvec bnew(9);
  arma::vec beta0(10);
  arma::rowvec Ygamma(m);
  arma::rowvec gammanew(q);
  //arma::cube L(max(Ds),q,S);
  //arma::cube P(max(Ds),q,S);
  int lens, ns, nt;
  double fnow, fnew, likelinow, likelinew, xbeta;
  for(int i=0; i<S; i++){
    //cout<<i<<endl;
    lens=getind(E,i+1).n_cols;
    inds.set_size(1,lens);
    inds = getind(E,i+1);
    Gammai.set_size(Ds(i),q);
    Gammai = getGammai(Gamma,i,Ds(i),q);
    C.set_size(Ds(i), Ds(i));
    C = updateC(Gammai,theta,tau);
    if(Ds(i)>1){
      Cminus.set_size(Ds(i)-1, Ds(i)-1);
      b.set_size(1,Ds(i)-1);
      temp2.set_size(1,Ds(i));
      bnew.set_size(1,Ds(i)-1);
    }
    for(int j=0; j<Ds(i); j++){
      //cout<<j<<endl;
      ns = getind(R.row(i),j+1).n_cols;
      indi.set_size(1,ns);
      indi = getind(R.row(i),j+1);
      beta0.set_size(ns,1);
      temp.set_size(1,ns);
      for(int k = 0; k<ns; k++){
        //k in R is indi(k) now
        temp(k) =  findT(D.col(indi(k)).t());
      }
      nt = Nuni(temp);
      indT.set_size(1,nt);
      indT = unique(temp);
      if(Ds(i)==1){
        //need to change is we change the kernel function
        fnow = 1;
      }
      else{
        Cminus = removeii(C,j+1,j+1);
        b = removei(C.row(j),j+1);
        temp1 = b*inverse_via_chol(Cminus)*b.t();
        fnow = kernelC(Gammai.row(j),Gammai.row(j),theta,tau) - temp1(0,0);
      }
      likelinow = 0;
      for(int k = 0; k<lens;k++){
        //k in R is now inds(k)
        xbeta = sum(X.row(inds(k))%Beta.row(i));
        Ygamma = Y.row(inds(k)) -  xbeta;
        beta0 = rowsome(Z,indi)* Gammai.row(j).t();
        likelinow += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);

      }
      //update Gamma one element by one element
      for(int l=0; l<q; l++){
        gammanew = Gammai.row(j);
        gammanew(l) = gammanew(l) + step(j,l,i);
        likelinew = 0;
        for(int k=0; k<lens; k++){
          //k in R is now inds(k)
          xbeta = sum(X.row(inds(k))%Beta.row(i));
          Ygamma = Y.row(inds(k)) -  xbeta;
          beta0 = rowsome(Z,indi)* gammanew.t();
          likelinew += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);
          munew = putvalue(munew,inds(k),indi,
                           computemu(mu,inds(k),indi, rowsome(Z,indi), Gammai.row(j), gammanew));


        }
        if(Ds(i)==1){
          fnew = 1;
          //need to change is we change the kernel function
        }
        else{
          for(int k=0; k<Ds(i);k++){
            temp2(k) = kernelC(Gammai.row(k), gammanew,theta,tau);
          }
          bnew = removei(temp2, j+1);
          temp1 = bnew*inverse_via_chol(Cminus)*bnew.t();
          fnew = kernelC(gammanew,gammanew,theta, tau) - temp1(0,0);
        }
        //L(j,l,i) = likelinew - likelinow;
        //P(j,l,i) = log(fnew/fnow);
        if((lens==0) || (ns==0)){
          likelinew = 0;
          likelinow = 0;
        }
        if((log(fnew/fnow) + likelinew - likelinow)>log(runif(j,l,i))){
          if(Ds(i)>1){
            fnow = fnew;
            C.row(j) = temp2;
            C.col(j) = temp2.t();
            C(j,j) = kernelC(gammanew,gammanew,theta,tau);
          }
          mu = munew;
          likelinow = likelinew;
          Gammai.row(j) = gammanew;
        }
        else{
          munew = mu;
        }
      }
    }
    //Gamma[,,i] = Gammai
    Gamma = putGamma(Gamma, Gammai, i, q);
  }
  ans["mu"]= mu;
  ans["Gamma"]= Gamma;
  //ans["L"] = L;
  //ans["P"] = P;
  return(ans);
}


List updateE(arma::mat Beta,arma::cube Gamma,arma::rowvec w,
             arma::mat X, arma::mat Y,arma::mat Z, arma::mat delta,
             arma::rowvec E,arma::mat R,double S,arma::rowvec Ds,
             arma::mat mu,arma::mat mustar,double sigma,arma::rowvec c,
             int n, int m, int T0, int p, int q, arma::mat D){
  List ans;
  arma::rowvec pro(S);
  arma::vec beta0(m);
  arma::rowvec gammaij(q);
  arma::rowvec indj(m);
  arma::vec mustarnew(T0);
  int len;
  arma::mat Zj(m,q);
  double likeli, va;
  arma::vec temp3(m);
  for(int i=0; i<n;i++){
    //compute the likeli
    pro = w;
    for(int k=0; k<S; k++){
      beta0 = matrix(X.row(i),m)*Beta.row(k).t() ;

      for(int j=0; j<Ds(k); j++){
        len=getind(R.row(k),j+1).n_cols;
        indj.set_size(1,len);
        if(len == 0){
          continue;
        }
        indj = getind(R.row(k),j+1);
        Zj.set_size(len,q);
        Zj = rowsome(Z,indj);
        gammaij = getgamma(Gamma, k,j);
        temp3.set_size(len,1);
        temp3 = Zj*gammaij.t();
        beta0 = minusvalue(beta0.t(), 0, indj, -temp3.t()).t();
      }




      likeli = dmvnrm_arma(Y.row(i),beta0.t(),sigma, true);
      mustarnew = mustar.row(i).t() -cDmu(c(1),D*mu.row(i).t()) + cDmu(c(1), D*beta0);
      for(int t=0; t<T0; t++){
        va = min(5.0, mustarnew(t));
        va = max(va, -5.0);
        likeli += R::pnorm(va,0,1,true,true)*delta(i,t) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(i,t));
      }
      pro(k) = log(pro(k)) + likeli;
    }


    pro = pro -max(pro);
    pro = exp(pro);
    pro = pro/sum(pro);
    //assign indicator
    E(i) =  rmunoim(pro);
    //update mu and mustar
    beta0 = matrix(X.row(i),m)*Beta.row(E(i)-1).t() ;
    for(int j=0; j<Ds(E(i)-1); j++){
      len=getind(R.row(E(i)-1),j+1).n_cols;
      indj.set_size(1,len);
      indj = getind(R.row(E(i)-1),j+1);
      Zj.set_size(len,q);
      Zj = rowsome(Z,indj);
      gammaij = getgamma(Gamma, E(i)-1,j);
      temp3.set_size(len,1);
      temp3 = Zj*gammaij.t();
      beta0 = minusvalue(beta0.t(), 0, indj, -temp3.t()).t();
    }
    mustarnew = mustar.row(i).t() -cDmu(c(1),D*mu.row(i).t()) + cDmu(c(1), D*beta0);
    mu.row(i) = beta0.t();
    mustar.row(i) = mustarnew.t();
  }
  ans["E"] = E;
  ans["Ds"] = Ds;
  ans["mu"] = mu;
  ans["mustar"] = mustar;
  return(ans);
}


List updateEnomiss(arma::mat Beta, arma::cube Gamma, arma::rowvec w,
                   arma::mat X, arma::mat Y,arma::mat Z,
                   arma::mat delta,arma::rowvec E,arma::mat R,
                   double S,arma::rowvec Ds,arma::mat mu,double sigma,
                   int n, int m, int T0, int p, int q, arma::mat D){
  List ans;
  arma::rowvec pro(S);
  arma::vec beta0(m);
  arma::rowvec gammaij(q);
  arma::rowvec indj(m);
  int len;
  arma::mat Zj(m,q);
  double likeli;
  arma::vec temp3(m);
  for(int i=0; i<n;i++){
    //compute the likeli
    pro = w;
    for(int k=0; k<S; k++){
      beta0 = matrix(X.row(i),m)*Beta.row(k).t() ;

      for(int j=0; j<Ds(k); j++){
        len=getind(R.row(k),j+1).n_cols;
        indj.set_size(1,len);
        if(len == 0){
          continue;
        }
        indj = getind(R.row(k),j+1);
        Zj.set_size(len,q);
        Zj = rowsome(Z,indj);
        gammaij = getgamma(Gamma, k,j);
        temp3.set_size(len,1);
        temp3 = Zj*gammaij.t();
        beta0 = minusvalue(beta0.t(), 0, indj, -temp3.t()).t();
      }

      likeli = dmvnrm_arma(Y.row(i),beta0.t(),sigma, true);

      pro(k) = log(pro(k)) + likeli;
    }


    pro = pro -max(pro);
    pro = exp(pro);
    pro = pro/sum(pro);
    //assign indicator
    E(i) =  rmunoim(pro);
    //update mu and mustar
    beta0 = matrix(X.row(i),m)*Beta.row(E(i)-1).t() ;
    for(int j=0; j<Ds(E(i)-1); j++){
      len=getind(R.row(E(i)-1),j+1).n_cols;
      indj.set_size(1,len);
      indj = getind(R.row(E(i)-1),j+1);
      Zj.set_size(len,q);
      Zj = rowsome(Z,indj);
      gammaij = getgamma(Gamma, E(i)-1,j);
      temp3.set_size(len,1);
      temp3 = Zj*gammaij.t();
      beta0 = minusvalue(beta0.t(), 0, indj, -temp3.t()).t();
    }
    mu.row(i) = beta0.t();
  }
  ans["E"] = E;
  ans["Ds"] = Ds;
  ans["mu"] = mu;
  return(ans);
}

arma::mat updatemu(arma::mat R,arma::mat Z,arma::mat X,
                   arma::cube Gamma,arma::rowvec K,arma::mat Beta,
                   arma::rowvec E,
                   int m, int n,int p, int q){
  arma::mat mu(n,m);
  arma::rowvec Xi(p);
  arma::rowvec Betai(p);
  arma::rowvec gammaij(q);
  arma::mat mui_temp(1,1);
  arma::mat Zj = Z;
  arma::rowvec indj(m);
  arma::vec temp1(m);
  int len;
  for(int i=0; i<n; i++){
    Xi = X.row(i);
    Betai = Beta.row(E(i)-1);
    mui_temp = sum(Xi%Betai);
    mu.row(i) = arma::repmat(mui_temp,1,m);
    for(int j=0; j<K(E(i)-1); j++){
      gammaij = Gamma.slice(E(i)-1).row(j);
      len = getind(R.row(E(i)-1),j+1).n_cols;
      indj.set_size(1,len);
      indj = getind(R.row(E(i)-1),j+1);
      Zj.set_size(len,q);
      Zj=rowsome(Z,indj);
      temp1.set_size(len,1);
      temp1 = Zj*gammaij.t();
      mu = minusvalue(mu, i, indj, -temp1.t());
    }
  }
  return (mu);
}

arma::rowvec putvalue_vec(arma::rowvec x, arma::rowvec indi, arma::colvec value){
  int n = indi.n_cols;
  for(int i=0; i<n; i++){
    x(indi(i)) = value(i);
  }
  return (x);
}

arma::mat rbind(arma::mat x, arma::rowvec y, arma::rowvec z){
  int R = x.n_rows+2;
  int C = x.n_cols;
  arma::mat ans(R,C);
  for(int i=0; i<R-2; i++){
    ans.row(i) = x.row(i);
  }
  ans.row(R-2) = y;
  ans.row(R-1) = z;
  return (ans);
}

arma::mat rbind0(arma::mat x, arma::rowvec y){
  int R = x.n_rows+1;
  int C = x.n_cols;
  arma::mat ans(R,C);
  for(int i=0; i<R-1; i++){
    ans.row(i) = x.row(i);
  }
  ans.row(R-1) = y;
  return (ans);
}

arma::rowvec connect(arma::rowvec x, double y, double z){
  int n = x.n_cols;
  x.reshape(1,n+2);
  x(n) = y;
  x(n+1) = z;
  return(x);
}

arma::rowvec connect0(arma::rowvec x, double y){
  int n = x.n_cols;
  x.reshape(1,n+1);
  x(n) = y;
  return(x);
}

double det_via_chol(arma::mat A){
  int size=A.n_cols;
  arma::mat chol = choleskyDecompose(A);
  double ans = 1.0;

  for(int i=0; i<size; i++){
    ans *= chol(i,i);
  }
  return ans*ans;
}


List Split(arma::rowvec w,int K,arma::mat Gamma,arma::rowvec Beta,
           arma::mat X,arma::mat Y, arma::mat Z,arma::rowvec R,
           arma::mat delta,arma::rowvec mu, arma::mat mu_star,
           double c,double sigma_square, arma::mat C,
           double theta, double tau,
           int m, int n, int q, int T0, double hyper_delta = 1){
  List ans;
  int Knew = K + 1;
  int l = sampleint(K);
  double w_1 = w(l-1);
  int nid, Tl;
  double va;
  arma::mat temp(1,1);
  //this temp is used for repmat
  arma::rowvec gamma_1(q);
  nid=getind(R,l).n_cols;
  arma::rowvec l_label(nid);
  arma::rowvec indT(nid);
  l_label = getind(R,l);
  for(int i=0; i<nid; i++){
    indT(i) = ceil((l_label(i) + 1)/6.0) - 1;
  }
  int nt = Nuni(indT);
  arma::rowvec indT_unique(nt);
  indT_unique = unique(indT);
  arma::mat Y_use(n,nid);
  arma::mat Z_use(nid,m);
  Z_use = rowsome(Z,l_label);
  Y_use = colsome(Y,l_label);

  //re-arrange labels
  arma::rowvec R_new(m);

  for(int i=0; i<m; i++){
    if(R(i)>l){
      R_new(i) = R(i)-1;
    }
    else{
      R_new(i) = R(i);
    }
  }
  //compute the split transformation using moment matching principle
  arma::rowvec beta(q);
  for(int i=0; i<q; i++){
    beta(i) = R::rbeta(2,2);
  }
  double alpha = R::runif(0,1);
  double w_new_1 = w_1 * alpha;
  double w_new_2 = w_1 * (1-alpha);
  arma::rowvec pro(2);
  pro(0) = log(w_new_1);
  pro(1) = log(w_new_2);
  arma::rowvec gamma_new_1(q);
  arma::rowvec gamma_new_2(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    gamma_new_1(i) = gamma_1(i) - std::sqrt((1-alpha)/alpha) * beta(i);
    gamma_new_2(i) = gamma_1(i) + std::sqrt(alpha/(1-alpha)) * beta(i);
    ratio_propose += R::dbeta(beta(i),2,2,true);
  }

  //compute the proposal ratio
  arma::rowvec split_to_two(nid);
  split_to_two.fill(0);
  int n1 = 0;
  int n2 = 0;
  double ratio_like = 0;
  if(nid>0){
    for(int i=0; i<nid; i++){
      temp(0,0) = sum(gamma_new_1 % Z_use.row(i));
      pro(0) += dmvnrm_arma(Y_use.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_new_2 % Z_use.row(i));
      pro(1) += dmvnrm_arma(Y_use.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      if(R::runif(0,1)<pro(0)){
        split_to_two(i) = 1;
        ratio_propose += log(pro(0));
        n1 += 1;
        R_new(l_label(i)) = K;
      }
      else{
        ratio_propose += log(pro(1));
        n2 += 1;
        R_new(l_label(i)) = Knew;
      }

    }

    arma::rowvec ind1(n1);
    arma::rowvec ind2(n2);
    ind1 = getind(split_to_two, 1);
    ind2 = getind(split_to_two, 0);
    arma::rowvec mu_new(nid);
    arma::rowvec mu_now(nid);
    mu_now = colsome(mu, l_label);
    mu_new = mu_now;
    double likeli_now = 0;
    double likeli_new = 0;
    arma::mat Z_1(n1,q);
    arma::mat Z_2(n2,q);
    Z_1 = rowsome(Z_use, ind1);
    Z_2 = rowsome(Z_use, ind2);
    mu_new = putvalue_vec(mu_new, ind1, Z_1*gamma_new_1.t());
    mu_new = putvalue_vec(mu_new, ind2, Z_2*gamma_new_2.t());
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y_use.row(i), mu_now,  sigma_square, true);
      likeli_new += dmvnrm_arma(Y_use.row(i), mu_new,  sigma_square, true);
    }

    arma::mat mu_star_now(n,T0);
    arma::mat mu_star_new(n,T0);
    mu_star_new = mu_star;
    mu_star_now = mu_star;

    for(int j=0; j<n; j++){
      for(int i=0; i<nid; i++){
        mu_star_new(j,indT(i)) = mu_star_new(j,indT(i)) - c * mu_now(i) * (1.0/6) + c * mu_new(i) * (1.0/6);
      }
    }

    for(int j=0; j<n; j++){
      for(int i=0; i<nt; i++){
        Tl = indT_unique(i);
        va = min(5.0, mu_star_now(j,Tl) );
        va = max(va, -5.0);
        likeli_now += R::pnorm(va,0,1,true,true)*delta(j,Tl) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(j,Tl));
        va = min(5.0, mu_star_new(j,Tl) );
        va = max(va, -5.0);
        likeli_new += R::pnorm(va,0,1,true,true)*delta(j,Tl) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(j,Tl));
      }
    }

    ratio_like = likeli_new - likeli_now;
  }

  //compute the log Jacobian
  double log_Jacobian = log(w_1) - (q/2.0) * log(alpha*(1-alpha));

  //compute the likelihood ratio

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w_new_1) + (hyper_delta - 1 + n2) * log(w_new_2) -
    (hyper_delta - 1 + nid) * log(w_1) - log(R::beta(hyper_delta, K * hyper_delta));

  arma::mat Gamma_new(Knew,q);
  arma::mat C_new(Knew,Knew);
  arma::rowvec w_new(Knew);
  if(n1==0){
    Knew -= 1;
    Gamma_new.set_size(Knew,q);
    C_new.set_size(Knew,Knew);
    w_new.set_size(1,Knew);
    w_new = connect0(removei(w,l),w_new_2);
    Gamma_new = rbind0(removeiii(Gamma, l),gamma_new_2);
    for(int i=0; i<nid; i++){
      R_new(l_label(i)) -= 1;
    }

  }else if(n2==0){
    Knew -= 1;
    Gamma_new.set_size(Knew,q);
    C_new.set_size(Knew,Knew);
    w_new.set_size(1,Knew);
    w_new = connect0(removei(w,l),w_new_1);
    Gamma_new = rbind0(removeiii(Gamma, l),gamma_new_1);

  }else{
    w_new = connect(removei(w,l), w_new_1, w_new_2);
    Gamma_new = rbind(removeiii(Gamma, l),gamma_new_1, gamma_new_2);
  }
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C_new)) - log(det_via_chol(C));
  double ratio_prior = ratio_prior_gamma + ratio_prior_w;
  double birth_ratio = ratio_like + ratio_prior - log(Knew) + log_Jacobian - log(Knew) - ratio_propose;
  if(Knew == 2){
    birth_ratio -= log(2);
  }
  if(log(R::runif(0,1)) < birth_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = w_new;
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}


List Splitnomiss(arma::rowvec w,int K,
                 arma::mat Gamma,arma::rowvec Beta,
                 arma::mat X,arma::mat Y,arma::mat Z,arma::rowvec R,
                 arma::mat delta,arma::rowvec mu,double sigma_square,arma::mat C,
                 double theta, double tau,
                 int m, int n, int q, int T0, double hyper_delta = 1){
  List ans;
  int Knew = K + 1;
  int l = sampleint(K);
  double w_1 = w(l-1);
  int nid;
  arma::mat temp(1,1);
  //this temp is used for repmat
  arma::rowvec gamma_1(q);
  nid=getind(R,l).n_cols;
  arma::rowvec l_label(nid);
  arma::rowvec indT(nid);
  l_label = getind(R,l);
  for(int i=0; i<nid; i++){
    indT(i) = ceil((l_label(i) + 1)/6.0) - 1;
  }
  int nt = Nuni(indT);
  arma::rowvec indT_unique(nt);
  indT_unique = unique(indT);
  arma::mat Y_use(n,nid);
  arma::mat Z_use(nid,m);
  Z_use = rowsome(Z,l_label);
  Y_use = colsome(Y,l_label);

  //re-arrange labels
  arma::rowvec R_new(m);

  for(int i=0; i<m; i++){
    if(R(i)>l){
      R_new(i) = R(i)-1;
    }
    else{
      R_new(i) = R(i);
    }
  }
  //compute the split transformation using moment matching principle
  arma::rowvec beta(q);
  for(int i=0; i<q; i++){
    beta(i) = R::rbeta(2,2);
  }
  double alpha = R::runif(0,1);
  double w_new_1 = w_1 * alpha;
  double w_new_2 = w_1 * (1-alpha);
  arma::rowvec pro(2);
  pro(0) = log(w_new_1);
  pro(1) = log(w_new_2);
  arma::rowvec gamma_new_1(q);
  arma::rowvec gamma_new_2(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    gamma_new_1(i) = gamma_1(i) - std::sqrt((1-alpha)/alpha) * beta(i);
    gamma_new_2(i) = gamma_1(i) + std::sqrt(alpha/(1-alpha)) * beta(i);
    ratio_propose += R::dbeta(beta(i),2,2,true);
  }

  //compute the proposal ratio
  arma::rowvec split_to_two(nid);
  split_to_two.fill(0);
  int n1 = 0;
  int n2 = 0;
  double ratio_like = 0;
  if(nid>0){
    for(int i=0; i<nid; i++){
      temp(0,0) = sum(gamma_new_1 % Z_use.row(i));
      pro(0) += dmvnrm_arma(Y_use.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_new_2 % Z_use.row(i));
      pro(1) += dmvnrm_arma(Y_use.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      if(R::runif(0,1)<pro(0)){
        split_to_two(i) = 1;
        ratio_propose += log(pro(0));
        n1 += 1;
        R_new(l_label(i)) = K;
      }
      else{
        ratio_propose += log(pro(1));
        n2 += 1;
        R_new(l_label(i)) = Knew;
      }

    }

    arma::rowvec ind1(n1);
    arma::rowvec ind2(n2);
    ind1 = getind(split_to_two, 1);
    ind2 = getind(split_to_two, 0);
    arma::rowvec mu_new(nid);
    arma::rowvec mu_now(nid);
    mu_now = colsome(mu, l_label);
    mu_new = mu_now;
    double likeli_now = 0;
    double likeli_new = 0;
    arma::mat Z_1(n1,q);
    arma::mat Z_2(n2,q);
    Z_1 = rowsome(Z_use, ind1);
    Z_2 = rowsome(Z_use, ind2);
    mu_new = putvalue_vec(mu_new, ind1, Z_1*gamma_new_1.t());
    mu_new = putvalue_vec(mu_new, ind2, Z_2*gamma_new_2.t());
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y_use.row(i), mu_now,  sigma_square, true);
      likeli_new += dmvnrm_arma(Y_use.row(i), mu_new,  sigma_square, true);
    }



    ratio_like = likeli_new - likeli_now;
  }

  //compute the log Jacobian
  double log_Jacobian = log(w_1) - (q/2.0) * log(alpha*(1-alpha));

  //compute the likelihood ratio

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w_new_1) + (hyper_delta - 1 + n2) * log(w_new_2) -
    (hyper_delta - 1 + nid) * log(w_1) - log(R::beta(hyper_delta, K * hyper_delta));

  arma::mat Gamma_new(Knew,q);
  arma::mat C_new(Knew,Knew);
  arma::rowvec w_new(Knew);
  if(n1==0){
    Knew -= 1;
    Gamma_new.set_size(Knew,q);
    C_new.set_size(Knew,Knew);
    w_new.set_size(1,Knew);
    w_new = connect0(removei(w,l),w_new_2);
    Gamma_new = rbind0(removeiii(Gamma, l),gamma_new_2);
    for(int i=0; i<nid; i++){
      R_new(l_label(i)) -= 1;
    }

  }else if(n2==0){
    Knew -= 1;
    Gamma_new.set_size(Knew,q);
    C_new.set_size(Knew,Knew);
    w_new.set_size(1,Knew);
    w_new = connect0(removei(w,l),w_new_1);
    Gamma_new = rbind0(removeiii(Gamma, l),gamma_new_1);

  }else{
    w_new = connect(removei(w,l), w_new_1, w_new_2);
    Gamma_new = rbind(removeiii(Gamma, l),gamma_new_1, gamma_new_2);
  }
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C_new)) - log(det_via_chol(C));
  double ratio_prior = ratio_prior_gamma + ratio_prior_w;
  double birth_ratio = ratio_like + ratio_prior - log(Knew) + log_Jacobian - log(Knew) - ratio_propose;
  if(Knew == 2){
    birth_ratio -= log(2);
  }
  if(log(R::runif(0,1)) < birth_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = w_new;
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}



arma::rowvec sampleint2(int K){
  arma::rowvec ans(2);
  int a = sampleint(K);
  int b = sampleint(K-1);
  if(b>=a){
    b += 1;
    ans(0) = a;
    ans(1) =b;
  }
  else{
    ans(0)=b;
    ans(1) = a;
  }
  return (ans);
}


List Merge(arma::rowvec w,
           int K,
           arma::mat Gamma,
           arma::rowvec Beta,
           arma::mat X,
           arma::mat Y,
           arma::mat Z,
           arma::rowvec R,
           arma::mat delta,
           arma::rowvec mu,
           arma::mat mu_star,
           double c,
           double sigma_square,
           arma::mat C,
           double theta, double tau,
           int m, int n, int q, int T0, double hyper_delta = 1){
  List ans;
  int Knew = K-1;

  //Uniformly pick two labels that are potentially to be merged
  arma::rowvec id(2);
  id = sampleint2(K);
  int id1 = id(0);
  int id2 = id(1);
  double w1 = w(id1-1);
  double w2 = w(id2-1);
  arma::rowvec pro_init(2);
  arma::rowvec pro(2);
  pro_init(0) = log(w1);
  pro_init(1) = log(w2);
  arma::rowvec gamma_1(q);
  arma::rowvec gamma_2(q);
  arma::rowvec gamma_new(q);
  gamma_1 = Gamma.row(id1-1);
  gamma_2 = Gamma.row(id2-1);
  int Tl;
  double va;
  // Find those groups that are associated with the kth mixture
  int n1=getind(R,id1).n_cols;
  int n2=getind(R,id2).n_cols;
  arma::mat Z1(n,q);
  arma::mat Y1(n,n);
  arma::mat Z2(n,q);
  arma::mat Y2(n,n);
  arma::rowvec ind1(m);
  arma::rowvec ind2(m);
  int nid = n1+n2;
  if(n1>0){
    ind1.set_size(n1);
    ind1 = getind(R,id1);
    Z1.set_size(n1,q);
    Y1.set_size(n,n1);
    Z1 = rowsome( Z, ind1);
    Y1 = colsome(Y, ind1);
  }
  if(n2>0){
    ind2.set_size(n2);
    ind2 = getind(R,id2);
    Z2.set_size(n2,q);
    Y2.set_size(n,n2);
    Z2 = rowsome(Z,ind2);
    Y2 = colsome(Y, ind2);
  }

  //re-arrange labels
  arma::rowvec R_new(m);
  R_new = R;
  for(int i=0; i<m; i++){
    if((R(i)>id1) & (R(i)<id2)){
      R_new(i)= R(i) -1;
    }
    if((R(i)>id2)){
      R_new(i)= R(i) -2;
    }
    if((R(i)==id1) || (R(i)==id2)){
      R_new(i) = Knew;
    }
  }

  //Compute the merge transformation using moment-matching principle
  double w1_new = w1 + w2;
  arma::rowvec w_new(Knew);
  w_new(Knew - 1) = w1_new;
  for(int i=0; i<Knew-1; i++){
    if(i<id1-1){
      w_new(i) = w(i);
    }
    if((i>=id1-1) & (i<id2-1)){
      w_new(i) = w(i+1);
    }
    if(i>=id2-1){
      w_new(i) = w(i+2);
    }
  }

  double alpha = w1/w1_new;
  arma::mat Gamma_new(Knew,q);
  gamma_new = (w1 * gamma_1 + w2 * gamma_2)/w1_new;
  Gamma_new = rbind0(removeiii(removeiii(Gamma,id1), id2-1),gamma_new);

  arma::rowvec beta(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    beta(i) = (gamma_new(i) - gamma_1(i)) * std::sqrt(alpha/(1-alpha));
    if((beta(i)<1) & (beta(i)>0)){
      ratio_propose += R::dbeta(beta(i),2,2,true);
    }
  }
  //compute the proposal ratio
  double diff = 0;
  arma::mat temp(1,1);
  if(n1>0){
    for(int i=0; i<n1; i++){
      pro = pro_init;
      temp(0,0) = sum(gamma_1 % Z1.row(i));
      pro(0) += dmvnrm_arma(Y1.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_2 % Z1.row(i));
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      diff = pro(1) - pro(0);
      pro(0) += 0.00001*diff;
      pro(1) -= 0.00001*diff;
      ratio_propose += log(pro(0));
    }
  }
  if(n2>0){
    for(int i=0; i<n2; i++){
      pro = pro_init;
      temp(0,0) = sum(gamma_1 % Z2.row(i));
      pro(0) += dmvnrm_arma(Y2.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_2 % Z2.row(i));
      pro(1) += dmvnrm_arma(Y2.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      diff = pro(1) - pro(0);
      pro(0) += 0.00001*diff;
      pro(1) -= 0.00001*diff;
      ratio_propose += log(pro(1));
    }
  }

  //compute the log Jacobian of moment-matching transformation
  double log_Jacobian = log(w1_new) - (q/2.0) * log(alpha*(1-alpha));
  //compute the likelihood ratio
  double likeli_now = 0;
  double likeli_new = 0;
  arma::mat mu_star_now(n,T0);
  arma::mat mu_star_new(n,T0);
  mu_star_new = mu_star;
  mu_star_now = mu_star;
  if(n1>0){
    arma::rowvec mu_new1(n1);
    arma::rowvec mu_now1(n1);
    arma::rowvec indT1(n1);
    mu_now1 = colsome(mu,ind1);
    mu_new1 = gamma_new * Z1.t();
    for(int i=0; i<n1; i++){
      indT1(i) = ceil((ind1(i) + 1)/6.0) - 1;
    }
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y1.row(i) ,mu_now1,  sigma_square,true);
      likeli_new += dmvnrm_arma(Y1.row(i) ,mu_new1,  sigma_square,true);
    }
    for(int j=0; j<n; j++){
      for(int i=0; i<n1; i++){
        mu_star_new(j,indT1(i)) = mu_star_new(j,indT1(i)) - c * mu_now1(i) * (1.0/6) + c * mu_new1(i) * (1.0/6);
      }
    }

  }
  if(n2>0){
    arma::rowvec mu_new2(n2);
    arma::rowvec mu_now2(n2);
    arma::rowvec indT2(n2);
    mu_now2 = colsome(mu,ind2);
    mu_new2 = gamma_new * Z2.t();
    for(int i=0; i<n2; i++){
      indT2(i) = ceil((ind2(i) + 1)/6.0) - 1;
    }
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y2.row(i) ,mu_now2, sigma_square,true);
      likeli_new += dmvnrm_arma(Y2.row(i) ,mu_new2, sigma_square,true);
    }
    for(int j=0; j<n; j++){
      for(int i=0; i<n2; i++){
        mu_star_new(j,indT2(i)) = mu_star_new(j,indT2(i)) - c * mu_now2(i) * (1.0/6) + c * mu_new2(i) * (1.0/6);
      }
    }
  }
  int flag = 0;
  arma::rowvec indT_temp(m);
  if(nid>0){
    //compute indT_unique
    for(int i=0; i<m; i++){
      if((R(i)==id1) || (R(i)==id2)){
        indT_temp(flag) = ceil((i+1)/6.0) - 1;
        flag += 1;
      }
    }
    indT_temp.reshape(1,flag);
    int nt = Nuni(indT_temp);
    arma::rowvec indT_unique(nt);
    indT_unique = unique(indT_temp);
    for(int j=0; j<n; j++){
      for(int i=0; i<nt; i++){
        Tl = indT_unique(i);
        va = min(5.0, mu_star_now(j,Tl) );
        va = max(va, -5.0);
        likeli_now += R::pnorm(va,0,1,true,true)*delta(j,Tl) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(j,Tl));
        va = min(5.0, mu_star_new(j,Tl) );
        va = max(va, -5.0);
        likeli_new += R::pnorm(va,0,1,true,true)*delta(j,Tl) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(j,Tl));
      }
    }
  }
  double ratio_like = likeli_new - likeli_now;

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w1) + (hyper_delta - 1 + n2) * log(w2) -
    (hyper_delta - 1 + nid) * log(w1_new) - log(R::beta(hyper_delta, Knew * hyper_delta));
  arma::mat C_new(Knew, Knew);
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C)) - log(det_via_chol(C_new));
  double ratio_prior = ratio_prior_w + ratio_prior_gamma;
  double death_ratio = log(K)*2 + ratio_propose + ratio_like - ratio_prior - log_Jacobian;
  if(Knew ==1){
    death_ratio = death_ratio + log(2);
  }

  if(log(R::runif(0,1)) < death_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = w_new;
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}


List Mergenomiss(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec Beta,
                 arma::mat X,
                 arma::mat Y,
                 arma::mat Z,
                 arma::rowvec R,
                 arma::mat delta,
                 arma::rowvec mu,
                 double sigma_square,
                 arma::mat C,
                 double theta, double tau,
                 int m, int n, int q, int T0, double hyper_delta = 1){
  List ans;
  int Knew = K-1;

  //Uniformly pick two labels that are potentially to be merged
  arma::rowvec id(2);
  id = sampleint2(K);
  int id1 = id(0);
  int id2 = id(1);
  double w1 = w(id1-1);
  double w2 = w(id2-1);
  arma::rowvec pro_init(2);
  arma::rowvec pro(2);
  pro_init(0) = log(w1);
  pro_init(1) = log(w2);
  arma::rowvec gamma_1(q);
  arma::rowvec gamma_2(q);
  arma::rowvec gamma_new(q);
  gamma_1 = Gamma.row(id1-1);
  gamma_2 = Gamma.row(id2-1);
  // Find those groups that are associated with the kth mixture
  int n1=getind(R,id1).n_cols;
  int n2=getind(R,id2).n_cols;
  arma::mat Z1(n,q);
  arma::mat Y1(n,n);
  arma::mat Z2(n,q);
  arma::mat Y2(n,n);
  arma::rowvec ind1(m);
  arma::rowvec ind2(m);
  int nid = n1+n2;
  if(n1>0){
    ind1.set_size(n1);
    ind1 = getind(R,id1);
    Z1.set_size(n1,q);
    Y1.set_size(n,n1);
    Z1 = rowsome( Z, ind1);
    Y1 = colsome(Y, ind1);
  }
  if(n2>0){
    ind2.set_size(n2);
    ind2 = getind(R,id2);
    Z2.set_size(n2,q);
    Y2.set_size(n,n2);
    Z2 = rowsome(Z,ind2);
    Y2 = colsome(Y, ind2);
  }

  //re-arrange labels
  arma::rowvec R_new(m);
  R_new = R;
  for(int i=0; i<m; i++){
    if((R(i)>id1) & (R(i)<id2)){
      R_new(i)= R(i) -1;
    }
    if((R(i)>id2)){
      R_new(i)= R(i) -2;
    }
    if((R(i)==id1) || (R(i)==id2)){
      R_new(i) = Knew;
    }
  }

  //Compute the merge transformation using moment-matching principle
  double w1_new = w1 + w2;
  arma::rowvec w_new(Knew);
  w_new(Knew - 1) = w1_new;
  for(int i=0; i<Knew-1; i++){
    if(i<id1-1){
      w_new(i) = w(i);
    }
    if((i>=id1-1) & (i<id2-1)){
      w_new(i) = w(i+1);
    }
    if(i>=id2-1){
      w_new(i) = w(i+2);
    }
  }

  double alpha = w1/w1_new;
  arma::mat Gamma_new(Knew,q);
  gamma_new = (w1 * gamma_1 + w2 * gamma_2)/w1_new;
  Gamma_new = rbind0(removeiii(removeiii(Gamma,id1), id2-1),gamma_new);

  arma::rowvec beta(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    beta(i) = (gamma_new(i) - gamma_1(i)) * std::sqrt(alpha/(1-alpha));
    if((beta(i)<1) & (beta(i)>0)){
      ratio_propose += R::dbeta(beta(i),2,2,true);
    }
  }
  //compute the proposal ratio
  double diff = 0;
  arma::mat temp(1,1);
  if(n1>0){
    for(int i=0; i<n1; i++){
      pro = pro_init;
      temp(0,0) = sum(gamma_1 % Z1.row(i));
      pro(0) += dmvnrm_arma(Y1.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_2 % Z1.row(i));
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      diff = pro(1) - pro(0);
      pro(0) += 0.00001*diff;
      pro(1) -= 0.00001*diff;
      ratio_propose += log(pro(0));
    }
  }
  if(n2>0){
    for(int i=0; i<n2; i++){
      pro = pro_init;
      temp(0,0) = sum(gamma_1 % Z2.row(i));
      pro(0) += dmvnrm_arma(Y2.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      temp(0,0) = sum(gamma_2 % Z2.row(i));
      pro(1) += dmvnrm_arma(Y2.col(i).t() ,repmat(temp,1,n), sigma_square,true);
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      diff = pro(1) - pro(0);
      pro(0) += 0.00001*diff;
      pro(1) -= 0.00001*diff;
      ratio_propose += log(pro(1));
    }
  }

  //compute the log Jacobian of moment-matching transformation
  double log_Jacobian = log(w1_new) - (q/2.0) * log(alpha*(1-alpha));
  //compute the likelihood ratio
  double likeli_now = 0;
  double likeli_new = 0;

  if(n1>0){
    arma::rowvec mu_new1(n1);
    arma::rowvec mu_now1(n1);
    arma::rowvec indT1(n1);
    mu_now1 = colsome(mu,ind1);
    mu_new1 = gamma_new * Z1.t();
    for(int i=0; i<n1; i++){
      indT1(i) = ceil((ind1(i) + 1)/6.0) - 1;
    }
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y1.row(i) ,mu_now1,  sigma_square,true);
      likeli_new += dmvnrm_arma(Y1.row(i) ,mu_new1,  sigma_square,true);
    }


  }
  if(n2>0){
    arma::rowvec mu_new2(n2);
    arma::rowvec mu_now2(n2);
    arma::rowvec indT2(n2);
    mu_now2 = colsome(mu,ind2);
    mu_new2 = gamma_new * Z2.t();
    for(int i=0; i<n2; i++){
      indT2(i) = ceil((ind2(i) + 1)/6.0) - 1;
    }
    for(int i=0; i<n; i++){
      likeli_now += dmvnrm_arma(Y2.row(i) ,mu_now2, sigma_square,true);
      likeli_new += dmvnrm_arma(Y2.row(i) ,mu_new2, sigma_square,true);
    }
  }

  double ratio_like = likeli_new - likeli_now;

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w1) + (hyper_delta - 1 + n2) * log(w2) -
    (hyper_delta - 1 + nid) * log(w1_new) - log(R::beta(hyper_delta, Knew * hyper_delta));
  arma::mat C_new(Knew, Knew);
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C)) - log(det_via_chol(C_new));
  double ratio_prior = ratio_prior_w + ratio_prior_gamma;
  double death_ratio = log(K)*2 + ratio_propose + ratio_like - ratio_prior - log_Jacobian;
  if(Knew ==1){
    death_ratio = death_ratio + log(2);
  }

  if(log(R::runif(0,1)) < death_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = w_new;
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}


List RJi(arma::rowvec w,
         int K,
         arma::mat Gamma,
         arma::rowvec Beta,
         arma::mat X,
         arma::mat Y,
         arma::mat Z,
         arma::rowvec R,
         arma::mat delta,
         arma::rowvec mu,
         arma::mat mu_star,
         double c,
         double sigma_square,
         arma::mat C,
         double theta, double tau,
         int m, int n, int q, int T0, double hyper_delta = 1){
  arma::mat diff(n,m);
  diff = arma::repmat(X*Beta.t(),1,m);
  Y = Y - diff;
  mu = mu - diff.row(0);
  double merge_or_split = R::runif(0,1);
  if (((K == 1) || (merge_or_split < 0.5)) & (K<10)){
    return (
        Split(w, K, Gamma, Beta, X, Y, Z,
              R, delta, mu, mu_star, c, sigma_square,
              C, theta, tau, m, n, q, T0, hyper_delta)
    );
  }
  else{
    return (
        Merge(w, K, Gamma, Beta, X, Y, Z,
              R, delta, mu, mu_star, c, sigma_square,
              C, theta, tau, m, n, q, T0, hyper_delta)
    );
  }
}




List RJinomiss(arma::rowvec w,
               int K,
               arma::mat Gamma,
               arma::rowvec Beta,
               arma::mat X,
               arma::mat Y,
               arma::mat Z,
               arma::rowvec R,
               arma::mat delta,
               arma::rowvec mu,
               double sigma_square,
               arma::mat C,
               double theta, double tau,
               int m, int n, int q, int T0, double hyper_delta = 1){
  arma::mat diff(n,m);
  diff = arma::repmat(X*Beta.t(),1,m);
  Y = Y - diff;
  mu = mu - diff.row(0);
  double merge_or_split = R::runif(0,1);
  if (((K == 1) || (merge_or_split < 0.5)) & (K<10)){
    return (
        Splitnomiss(w, K, Gamma, Beta, X, Y, Z,
                    R, delta, mu,  sigma_square,
                    C, theta, tau, m, n, q, T0, hyper_delta)
    );
  }
  else{
    return (
        Mergenomiss(w, K, Gamma, Beta, X, Y, Z,
                    R, delta, mu, sigma_square,
                    C, theta, tau, m, n, q, T0, hyper_delta)
    );
  }
}



arma::rowvec updateRi(arma::rowvec w,
                      arma::mat Gamma,
                      arma::rowvec Beta,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::rowvec mu,
                      arma::mat mu_star,
                      double c,
                      double sigma_square,
                      int K, arma::mat X,
                      int m, int n, int q){
  arma::rowvec R(m);
  arma::mat diff(n,m);
  diff = arma::repmat(X*Beta.t(),1,m);
  Y = Y - diff;
  mu = mu - diff.row(0);
  double mu_s, va, likeli;
  int T_now;
  arma::mat mu_new(1,1);
  arma::rowvec pro(K);
  arma::colvec mu_star_s(n);
  arma::colvec mu_star_new_s(n);
  arma::rowvec Gammai(q);
  for(int s=0; s<m; s++){
    mu_s = mu(s);
    T_now = ceil((s + 1)/6.0) - 1;
    mu_star_s = mu_star.col(T_now);
    pro = log(w)/log_e;
    for(int k=0; k<K; k++){
      Gammai = Gamma.row(k);
      mu_new(0,0) = sum(Gammai % Z.row(s));
      mu_star_new_s = mu_star_s + c*(mu_new(0,0) - mu_s)*(1.0/6.0);
      likeli = 0;
      likeli += dmvnrm_arma(Y.col(s).t(), arma::repmat(mu_new,1,n), sigma_square, true);
      for(int i=0; i<n; i++){
        va = min(5.0, mu_star_new_s(i));
        va = max(va, -5.0);
        likeli += R::pnorm(va,0,1,true, true)* delta(i,T_now) + log(1-R::pnorm(va,0,1,true, false))*(1-delta(i,T_now));
      }
      pro(k) += likeli;
    }
    pro = pro -max(pro);
    pro = exp(pro);
    pro = pro/sum(pro);
    //assign indicator

    R(s) =  rmunoim(pro);
  }
  return (R);
}


arma::rowvec updateRinomiss(arma::rowvec w,
                            arma::mat Gamma,
                            arma::rowvec Beta,
                            arma::mat Y,
                            arma::mat Z,
                            arma::mat delta,
                            arma::rowvec mu,
                            double sigma_square,
                            int K, arma::mat X,
                            int m, int n, int q){
  arma::rowvec R(m);
  arma::mat diff(n,m);
  diff = arma::repmat(X*Beta.t(),1,m);
  Y = Y - diff;
  mu = mu - diff.row(0);
  double mu_s,  likeli;
  arma::mat mu_new(1,1);
  arma::rowvec pro(K);
  arma::rowvec Gammai(q);
  for(int s=0; s<m; s++){
    mu_s = mu(s);
    pro = log(w)/log_e;
    for(int k=0; k<K; k++){
      Gammai = Gamma.row(k);
      mu_new(0,0) = sum(Gammai % Z.row(s));
      likeli = 0;
      likeli += dmvnrm_arma(Y.col(s).t(), arma::repmat(mu_new,1,n), sigma_square, true);
      pro(k) += likeli;
    }
    pro = pro -max(pro);
    pro = exp(pro);
    pro = pro/sum(pro);
    //assign indicator

    R(s) =  rmunoim(pro);
  }
  return (R);
}


arma::rowvec sample_int(int size, arma::rowvec pro){
  arma::rowvec ans(size);
  for(int i=0; i<size; i++){
    ans(i) = rmunoim(pro);
  }
  return (ans);
}



arma::mat updateR(arma::mat w,
                  arma::cube Gamma,
                  arma::mat Beta,
                  arma::mat Y,
                  arma::mat Z,
                  arma::mat delta,
                  arma::mat mu,
                  arma::mat mu_star,
                  double c, int S,
                  double sigma_square,
                  arma::rowvec K,
                  arma::rowvec E,
                  arma::mat X,
                  int m, int n, int q, int p, int T0){
  arma::mat R(S,m);
  arma::rowvec wi(10);
  arma::mat Gammai(10,q);
  arma::rowvec Betai(p);
  arma::rowvec indi(n);
  arma::mat Yi(n,m);
  arma::mat deltai(n,T0);
  arma::rowvec mui(m);
  arma::mat mu_stari(n,T0);
  arma::rowvec Ri(m);
  arma::mat Xi(n,p);

  int len, Ki;
  for(int i=0; i<S; i++){
    wi.set_size(1,K(i));
    wi = colsome(w.row(i),ton(K(i)));
    len=getind(E,i+1).n_cols;
    if(len==0){
      R.row(i) = sample_int(m, wi);
    }
    else{
      Gammai.set_size(K(i),q);
      Gammai = getGammai(Gamma,i,K(i),q);
      Betai = Beta.row(i);
      indi.set_size(1,len);
      indi = getind(E,i+1);
      Yi.set_size(len,m);
      Yi = rowsome(Y,indi);
      deltai.set_size(len,T0);
      deltai = rowsome(delta,indi);
      mui = mu.row(indi(0));
      mu_stari.set_size(len,T0);
      mu_stari = rowsome(mu_star, indi);
      Ki = K(i);
      Xi.set_size(len,p);
      Xi = rowsome(X,indi);
      R.row(i) = updateRi(wi,Gammai, Betai, Yi, Z, deltai, mui, mu_stari, c, sigma_square, Ki,
            Xi, m, len, q);
    }


  }
  return (R);
}


arma::mat updateRnomiss(arma::mat w,
                        arma::cube Gamma,
                        arma::mat Beta,
                        arma::mat Y,
                        arma::mat Z,
                        arma::mat delta,
                        arma::mat mu,
                        int S,
                        double sigma_square,
                        arma::rowvec K,
                        arma::rowvec E,
                        arma::mat X,
                        int m, int n, int q, int p, int T0){
  arma::mat R(S,m);
  arma::rowvec wi(10);
  arma::mat Gammai(10,q);
  arma::rowvec Betai(p);
  arma::rowvec indi(n);
  arma::mat Yi(n,m);
  arma::mat deltai(n,T0);
  arma::rowvec mui(m);
  arma::rowvec Ri(m);
  arma::mat Xi(n,p);

  int len, Ki;
  for(int i=0; i<S; i++){
    wi.set_size(1,K(i));
    wi = colsome(w.row(i),ton(K(i)));
    len=getind(E,i+1).n_cols;
    if(len==0){
      R.row(i) = sample_int(m, wi);
    }
    else{
      Gammai.set_size(K(i),q);
      Gammai = getGammai(Gamma,i,K(i),q);
      Betai = Beta.row(i);
      indi.set_size(1,len);
      indi = getind(E,i+1);
      Yi.set_size(len,m);
      Yi = rowsome(Y,indi);
      deltai.set_size(len,T0);
      deltai = rowsome(delta,indi);
      mui = mu.row(indi(0));
      Ki = K(i);
      Xi.set_size(len,p);
      Xi = rowsome(X,indi);
      R.row(i) = updateRinomiss(wi,Gammai, Betai, Yi, Z, deltai, mui,  sigma_square, Ki,
            Xi, m, len, q);
    }


  }
  return (R);
}

arma::mat updatemustar(arma::mat mu,
                       arma::rowvec c,
                       int n, int T0, arma::mat D){
  arma::mat mustar(n,T0);
  arma::vec temp2(T0);
  for(int i=0; i<n;i++){
    temp2 = cDmu(c(1),D*mu.row(i).t()) + c(0);
    mustar.row(i) = temp2.t();
  }
  return (mustar);
}




double rtrun_norm( double mu, double sigma, double a, double b){
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;
  alpha_cdf = R::pnorm(alpha,0,1,true,false);
  beta_cdf = R::pnorm(beta,0,1,true,false);
  if(alpha_cdf>0.9999){
    return (R::runif(a,a+sigma));
  }else if(beta_cdf<0.0001){
    return (R::runif(b-sigma,b));
  }else{
    u = R::runif(0,1);
    xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
    xi = R::qnorm5(xi_cdf,0,1,true,false);
    x = mu + sigma * xi;

    return (x);
  }

}


arma::mat updateZstar(arma::mat mu_star, arma::mat delta,
                      int n, int T0){
  arma::mat Z(n,T0);
  for(int i=0; i<n; i++){
    for(int t=0; t<T0; t++){
      if(delta(i,t) == 0){
        Z(i,t) = rtrun_norm(mu_star(i,t),1,-Inf,0);
      }else{
        Z(i,t) = rtrun_norm(mu_star(i,t),1,0,Inf);
      }
    }
  }
  return (Z);
}


arma::rowvec rmvnorm_row(arma::rowvec mean, arma::mat var){
  int n = mean.n_cols;
  arma::rowvec ans(n);
  arma::rowvec Z(n);
  arma::mat R(n,n);
  arma::mat tmp(n,1);
  R = choleskyDecompose(var);
  for(int i=0; i<n; i++){
    Z(i) = R::rnorm(0,1);
  }
  tmp = R.t() * Z.t();
  ans = mean + tmp.t();
  //mu + R^TZ
  return (ans);

}



arma::rowvec updatec(arma::mat Zstar, arma::mat mu, arma::mat D,
                     double sigmac, double sigma_square, int n, int T0){
  double sumxx = 0;
  double sumxy = 0;
  double sumx = 0, sumy=0;
  arma::rowvec xi(T0);
  arma::rowvec Yc(T0);
  for(int i=0; i<n; i++){
    xi = (D * mu.row(i).t()).t();
    Yc = Zstar.row(i);
    sumx += sum(xi);
    sumxx += sum(xi % xi );
    sumxy += sum(xi % Yc);
    sumy += sum(Yc);
  }
  arma::mat Sig(2,2);
  arma::mat Sign(2,2);
  arma::rowvec exptation(2);
  exptation(0) = sumy/sigma_square;
  exptation(1) = sumxy/sigma_square;
  arma::rowvec ans(2);

  Sig(0,0) = n*T0;
  Sig(0,1) = sumx;
  Sig(1,0) = sumx;
  Sig(1,1) = sumxx;


  Sign = inverse_via_chol(Sig/sigma_square + eye(2)/sigmac);
  ans = rmvnorm_row(exptation*Sign, Sign);

  return (ans);

}


List updateBetanoDPP(arma::mat X,
                     arma::mat Y,
                     arma::mat Z,
                     arma::mat delta,
                     arma::mat Beta,
                     arma::cube Gamma,
                     arma::rowvec E,
                     arma::mat R,
                     double S,
                     arma::rowvec Ds,
                     arma::mat mustar,
                     arma::mat mu,
                     double sigma,
                     arma::rowvec c,
                     arma::mat step,
                     arma::mat runif,
                     int n, int m, int T0, int p, int q, arma::mat D){
  List ans;
  arma::rowvec Betanew(p);
  arma::vec mustarnewj(T0);
  arma::vec munewj(m);
  arma::mat mustarnew(n,T0);
  arma::mat munew(n,m);
  arma::mat Ybeta(n,m);
  arma::rowvec indj(m);
  arma::mat Zj(m,q);
  arma::rowvec gammaij(q);
  arma::rowvec betai(p);
  arma::rowvec indi(n);
  arma::vec beta0(m);
  int l;
  double likelinow, likelinew,va;
  munew = mu;
  mustarnew = mustar;
  Ybeta = Y;
  for(int i=0; i<n; i++){
    for(int j=1; j<=Ds(E(i)-1);j++){
      l=getind(R.row(E(i)-1),j).n_cols;
      indj.reshape(1,l);
      indj = getind(R.row(E(i)-1),j);
      Zj.reshape(l,q);
      Zj = rowsome(Z,indj);
      gammaij = Gamma.slice(E(i)-1).row(j-1);
      Ybeta = minusvalue(Ybeta,i,indj, (Zj*gammaij.t()).t());
    }
  }
  //now update beta one by one row

  for(int i=0; i<S;i++){
    likelinow = 0;
    betai = Beta.row(i);
    l=getind(E,i+1).n_cols;
    indi.reshape(1,l);
    indi = getind(E,i+1);
    for(int j=0; j<l; j++){
      //j in R is indi(j) now
      beta0 = matrix(X.row(indi(j)),m) * betai.t();
      likelinow += dmvnrm_arma(Ybeta.row(indi(j)),  beta0.t(),   sigma , true);
      for(int t=0; t<T0; t++){
        va = min(5.0, mustar(indi(j),t));
        va = max(va, -5.0);
        likelinow += R::pnorm(va,0,1,true, true)* delta(indi(j),t) + log(1-R::pnorm(va,0,1,true, false))*(1-delta(indi(j),t));
      }
    }
    //update Beta one by one

    for(int j=0; j<p;j++){
      Betanew = Beta.row(i);
      Betanew(j) = Betanew(j) + step(i,j);
      likelinew = 0;
      for(int k=0; k<l;k++){
        //k in R is now indi(k)
        beta0 = matrix(X.row(indi(k)),m)*Betanew.t();
        likelinew += dmvnrm_arma(Ybeta.row(indi(k)),  beta0.t(), sigma , true);
        munewj = mu.row(indi(k)).t() + matrix(X.row(indi(k)),m) * (Betanew.t()-Beta.row(i).t());
        mustarnewj = mustar.row(indi(k)).t() + cDmu(c(1),D*(munewj)) - cDmu(c(1),D*mu.row(indi(k)).t());
        mustarnew.row(indi(k)) = mustarnewj.t();
        munew.row(indi(k)) = munewj.t();
        for(int t=0; t<T0; t++){
          va = min(5.0, mustarnewj(t));
          va = max(va, -5.0);
          likelinew += R::pnorm(va,0,1,true,true)*delta(indi(k),t) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(indi(k),t));
        }
      }
      //compute prior
      if(l==0){
        likelinew = 0;
        likelinow = 0;
      }
      if((likelinew - likelinow + R::dnorm4(Betanew(j),0,100,true) - R::dnorm4(Beta(i,j),0,100,true))>log(runif(i,j))){
        likelinow = likelinew;
        Beta.row(i) = Betanew;
        mustar = mustarnew;
        mu = munew;
      }
    }
  }
  ans["Beta"] = Beta;
  ans["mu"]= mu;
  ans["mustar"] = mustar;
  return (ans);
}



List updateGammanoDPP(arma::mat X,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::mat Beta,
                      arma::cube Gamma,
                      arma::rowvec E,
                      arma::mat R,
                      double S,
                      arma::rowvec Ds,
                      arma::mat mu,
                      arma::mat mustar,
                      double sigma,
                      arma::rowvec c,
                      arma::cube step,
                      arma::cube runif,
                      int n, int m, int T0, int p, int q, arma::mat D){
  List ans;
  arma::mat mustarnew(n, T0);
  arma::mat munew(n, m);
  arma::rowvec inds(n);
  arma::rowvec indi(m);
  arma::rowvec indT(T0);
  arma::mat Gammai(max(Ds),q);
  mustarnew = mustar;
  arma::rowvec temp(T0);
  munew =mu;
  arma::mat temp1(1,1);
  arma::rowvec b(9);
  arma::rowvec temp2(10);
  arma::rowvec bnew(9);
  arma::vec beta0(10);
  arma::rowvec Ygamma(m);
  arma::mat Zi(m,q);
  arma::rowvec gammanew(q);
  //arma::cube L(max(Ds),q,S);
  //arma::cube P(max(Ds),q,S);
  int lens, ns, nt;
  double likelinow, likelinew, xbeta, va;
  for(int i=0; i<S; i++){
    //cout<<i<<endl;
    lens=getind(E,i+1).n_cols;
    inds.reshape(1,lens);
    inds = getind(E,i+1);
    Gammai.reshape(Ds(i),q);
    Gammai = getGammai(Gamma,i,Ds(i),q);
    for(int j=0; j<Ds(i); j++){
      //cout<<j<<endl;
      ns = getind(R.row(i),j+1).n_cols;
      indi.reshape(1,ns);
      indi = getind(R.row(i),j+1);
      beta0.reshape(ns,1);
      Zi.reshape(ns,q);
      temp.reshape(1,ns);
      for(int k = 0; k<ns; k++){
        //k in R is indi(k) now
        temp(k) =  findT(D.col(indi(k)).t());
      }
      nt = Nuni(temp);
      indT.reshape(1,nt);
      indT = unique(temp);
      likelinow = 0;
      for(int k = 0; k<lens;k++){
        //k in R is now inds(k)
        xbeta = sum(X.row(inds(k))%Beta.row(i));
        Ygamma = Y.row(inds(k)) -  xbeta;
        beta0 = rowsome(Z,indi)* Gammai.row(j).t();
        likelinow += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);
        for(int t=0; t<nt; t++){
          //t in R is now indT(t)
          va = min(5.0, mustar(inds(k),indT(t)));
          va = max(va, -5.0);
          likelinow += R::pnorm(va,0,1,true, true)*delta(inds(k),indT(t)) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(inds(k),indT(t)));
        }
      }
      //update Gamma one element by one element
      for(int l=0; l<q; l++){
        gammanew = Gammai.row(j);
        gammanew(l) = gammanew(l) + step(j,l,i);
        likelinew = 0;
        for(int k=0; k<lens; k++){
          //k in R is now inds(k)
          xbeta = sum(X.row(inds(k))%Beta.row(i));
          Ygamma = Y.row(inds(k)) -  xbeta;
          beta0 = rowsome(Z,indi)* gammanew.t();
          likelinew += dmvnrm_arma(colsome(Ygamma,indi), beta0.t(), sigma, true);
          munew = putvalue(munew,inds(k),indi,
                           computemu(mu,inds(k),indi, rowsome(Z,indi), Gammai.row(j), gammanew));
          mustarnew = putvalue(mustarnew,inds(k),indT,
                               computemustar(mustar, inds(k),indT, c, D, mu.row(inds(k)), munew.row(inds(k)), T0));
          for(int t=0; t<nt; t++){
            //t in R is now indT(t)
            va = min(5.0, mustarnew(inds(k),indT(t)));
            va = max(va, -5.0);
            likelinew += R::pnorm(va,0,1,true, true)*delta(inds(k),indT(t)) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(inds(k),indT(t)));

          }

        }
        if(lens == 0){
          likelinew =0;
          likelinow = 0;
        }
        if(ns == 0){
          likelinew = 0;
          likelinow = 0;
        }
        if((R::dnorm4(gammanew(l),0,100,true)- R::dnorm4(Gammai(j,l),0,100,true) + likelinew - likelinow)>log(runif(j,l,i))){
          mu = munew;
          mustar = mustarnew;
          likelinow = likelinew;
          Gammai.row(j) = gammanew;
        }
        else{
          munew = mu;
          mustarnew = mustar;
        }
      }
    }
    //Gamma[,,i] = Gammai
    Gamma = putGamma(Gamma, Gammai, i, q);
  }
  ans["mu"]= mu;
  ans["mustar"] = mustar;
  ans["Gamma"]= Gamma;
  //ans["L"] = L;
  //ans["P"] = P;
  return(ans);
}


arma::cube spcopyGamma2(arma::cube Gamma, arma::rowvec Ds, int S, int q){
  int max0;
  max0 = max(Ds);
  arma::cube ans(max0,q,S);
  ans.fill(0);
  for(int i=0; i<S; i++){
    for(int j=0; j<Ds(i); j++){
      ans.slice(i).row(j) = Gamma.slice(i).row(j);
    }
  }
  return (ans);
}

arma::cube spcopyGamma(arma::cube Gamma, int k, arma::mat Gammai, arma::rowvec Ds, int S, int q){
  //Gamma[,,i]<-Gammai
  //need to change dim
  //k starts from 0
  int max0;
  max0 = max(Ds);
  arma::cube ans(max0,q,S);
  ans.fill(0);
  for(int i=0; i<S; i++){
    if(i==k){
      for(int j=0; j<Ds(i); j++){
        ans.slice(i).row(j) = Gammai.row(j);
      }
    }
    else{
      for(int j=0; j<Ds(i); j++){
        ans.slice(i).row(j) = Gamma.slice(i).row(j);
      }
    }
  }
  return (ans);
}





List updateRnoDPP(arma::mat Beta,
                  arma::cube Gamma,
                  arma::mat X,
                  arma::mat Y,
                  arma::mat Z,
                  arma::mat delta,
                  arma::rowvec E,
                  arma::mat R,
                  int S,
                  arma::rowvec Ds,
                  arma::mat mu,
                  arma::mat mustar,
                  double sigma,
                  double c,
                  arma::rowvec Ms,
                  int n, int m, int T0, int p, int q, arma::mat D, int Dmax=10){
  List ans;
  arma::rowvec Ns(S);
  Ns.fill(0);
  for(int i=0;i<n;i++){
    Ns(E(i)-1)++;
  }
  arma::rowvec indi(n);
  arma::rowvec ns(m);
  arma::mat Yr = Y;
  double likeli, mu0now, mu0, mustarnew, va;
  int T1, max1, nnow;
  arma::rowvec pro(10);
  arma::vec temp1(m);
  arma::mat Gammai(10,q);
  arma::rowvec temp2(10);
  arma::mat temp3(max(Ds),q);
  arma::rowvec gammanew(q);
  arma::cube temp4=Gamma;
  for(int i=0; i<n; i++){
    temp1 = matrix(X.row(i),m)*Beta.row(E(i)-1).t();
    Yr.row(i) = Yr.row(i) - temp1.t();
  }
  for(int i=0; i<S; i++){
    Gammai.set_size(Ds(i),q);
    Gammai = getGammai(Gamma,i,Ds(i),q);
    indi.set_size(1,Ns(i));
    indi = getind(E,i+1);
    ns.set_size(1,Ds(i));
    pro.set_size(1,Ds(i)+1);
    for(int j=0; j<Ds(i); j++){
      ns(j) = sum(R.row(i)==(j+1));
    }
    for(int j=0; j<m; j++){
      //find the corresponding teeth
      T1 = findT(D.col(j).t());
      mu0now = sum(Z.row(j)%Gammai.row(R(i,j)-1));
      ns(R(i,j)-1)--;
      if(ns(R(i,j)-1)==0){
        Ds(i)--;
        if(R::runif(0,1)<((Ds(i)-1)/Ds(i))){
          ns(R(i,j)-1)++;
          Ds(i)++;
          continue;
        }
        temp2.set_size(1,Ds(i));
        temp2 = removei(ns, R(i,j));
        ns.set_size(1,Ds(i));
        pro.set_size(1,Ds(i)+1);
        ns = temp2;
        temp3.set_size(Ds(i),q);
        temp3 = removeiii(Gammai,R(i,j));
        Gammai.set_size(Ds(i),q);
        Gammai = temp3;
        for(int k=0; k<m; k++){
          if(R(i,j)<R(i,k)){
            R(i,k)--;
          }
        }
      }
      //compute prob
      pro= connect0(ns,Ms(i)/(Ds(i)+1));
      pro = log(pro);
      for(int k=0; k<Ds(i); k++){
        likeli = 0;
        mu0 = sum(Z.row(j)%Gammai.row(k));
        for(int l=0; l<Ns(i); l++){
          //l in R is now indi(l)
          if(Yr(indi(l),j)==Yr(indi(l),j)){
            likeli+= R::dnorm(Yr(indi(l),j),mu0, sigma, true);
          }

          mustarnew = mustar(indi(l), T1) - c*mu0now/6 + c*mu0/6;
          va = min(5.0, mustarnew);
          va = max(va, -5.0);
          likeli += R::pnorm(va,0,1,true,true)*delta(indi(l),T1) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(indi(l),T1));
        }
        pro(k) += likeli;
      }
      //for new cluster
      likeli = 0;
      for(int k=0; k<q; k++){
        gammanew(k) = R::rnorm(R::runif(-5,5), sqrt(5));
      }
      mu0=sum(Z.row(j)%gammanew);
      for(int l=0; l<Ns(i);l++){
        if(Yr(indi(l),j)==Yr(indi(l),j)){
          likeli+= R::dnorm(Yr(indi(l),j),mu0, sigma, true);
        }
        mustarnew = mustar(indi(l), T1) - mu0now/6 + mu0/6;
        va = min(5.0, mustarnew);
        va = max(va, -5.0);
        likeli += R::pnorm(va,0,1,true,true)*delta(indi(l),T1) + log(1-R::pnorm(va,0,1,true,false))*(1-delta(indi(l),T1));
      }
      pro(Ds(i)) += likeli;
      //compute the prior
      if(Ds(i)==Dmax){
        pro.reshape(1,Dmax);
      }
      pro = pro -max(pro);
      pro = exp(pro);
      pro = pro/sum(pro);
      //assign indicator
      R(i,j) =  rmunoim(pro);
      if(R(i,j)==(Ds(i)+1)){
        Ds(i)++;
        temp2.set_size(1,Ds(i));
        temp2 = connect0(ns,1);
        ns.set_size(1,Ds(i));
        ns = temp2;
        temp3.set_size(Ds(i),q);
        temp3 = rbind0(Gammai, gammanew);
        Gammai.set_size(Ds(i),q);
        Gammai = temp3;
      }
      else{
        ns(R(i,j)-1)++;
      }
    }
    //put Gammai back to Gamma
    nnow = Gamma.slice(0).n_rows;
    if(Ds(i)>nnow){
      temp4.set_size(Ds(i),q,S);
      temp4 = spcopyGamma(Gamma,i,Gammai,Ds,S,q);
      Gamma.set_size(Ds(i),q,S);
      Gamma = temp4;
    }
    else{
      Gamma = putGamma(Gamma, Gammai, i, q);
    }
  }
  max1 = max(Ds);
  nnow = Gamma.slice(0).n_rows;
  if(max1<nnow){
    temp4.set_size(max1,q,S);
    temp4 = spcopyGamma2(Gamma,Ds,S,q);
    Gamma.set_size(max1,q,S);
    Gamma = temp4;
  }
  ans["Gamma"] = Gamma;
  ans["Ds"] = Ds;
  ans["R"] = R;
  ans["Yr"] = Yr;
  return (ans);
}


List Split_empty(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec R,
                 arma::mat C,
                 double theta, double tau,
                 int m, int q, double hyper_delta = 1){
  List ans;
  int Knew = K + 1;
  int l = sampleint(K);
  double w_1 = w(l-1);
  int nid;
  arma::rowvec gamma_1(q);
  nid=getind(R,l).n_cols;
  arma::rowvec l_label(nid);
  l_label = getind(R,l);

  //re-arrange labels
  arma::rowvec R_new(m);

  for(int i=0; i<m; i++){
    if(R(i)>l){
      R_new(i) = R(i)-1;
    }
    else{
      R_new(i) = R(i);
    }
  }
  //compute the split transformation using moment matching principle
  arma::rowvec beta(q);
  for(int i=0; i<q; i++){
    beta(i) = R::rbeta(2,2);
  }
  double alpha = R::runif(0,1);
  double w_new_1 = w_1 * alpha;
  double w_new_2 = w_1 * (1-alpha);
  arma::rowvec pro(2);
  pro(0) = w_new_1;
  pro(1) = w_new_2;
  arma::rowvec gamma_new_1(q);
  arma::rowvec gamma_new_2(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    gamma_new_1(i) = gamma_1(i) - std::sqrt((1-alpha)/alpha) * beta(i);
    gamma_new_2(i) = gamma_1(i) + std::sqrt(alpha/(1-alpha)) * beta(i);
    ratio_propose += R::dbeta(beta(i),2,2,true);
  }

  //compute the proposal ratio
  int n1 = 0;
  int n2 = 0;
  if(nid>0){
    for(int i=0; i<nid; i++){
      if(R::runif(0,1)<pro(0)){
        ratio_propose += log(pro(0));
        n1 += 1;
        R_new(l_label(i)) = K;
      }
      else{
        ratio_propose += log(pro(1));
        n2 += 1;
        R_new(l_label(i)) = Knew;
      }
    }

  }

  //compute the log Jacobian
  double log_Jacobian = log(w_1) - (q/2.0) * log(alpha*(1-alpha));

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w_new_1) + (hyper_delta - 1 + n2) * log(w_new_2) -
    (hyper_delta - 1 + nid) * log(w_1) - log(R::beta(hyper_delta, K * hyper_delta));

  arma::mat Gamma_new(Knew,q);
  Gamma_new = rbind(removeiii(Gamma, l),gamma_new_1, gamma_new_2);
  arma::mat C_new(Knew,Knew);
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C_new)) - log(det_via_chol(C));
  double ratio_prior = ratio_prior_gamma + ratio_prior_w;
  double birth_ratio = ratio_prior - log(Knew) + log_Jacobian - log(Knew) - ratio_propose;
  if(Knew == 2){
    birth_ratio -= log(2);
  }
  if(log(R::runif(0,1)) < birth_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = connect(removei(w,l), w_new_1, w_new_2);
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}


List Merge_empty(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec R,
                 arma::mat C,
                 double theta, double tau,
                 int m,  int q, double hyper_delta = 1){
  List ans;
  int Knew = K-1;

  //Uniformly pick two labels that are potentially to be merged
  arma::rowvec id(2);
  id = sampleint2(K);
  int id1 = id(0);
  int id2 = id(1);
  double w1 = w(id1-1);
  double w2 = w(id2-1);
  arma::rowvec pro_init(2);
  arma::rowvec pro(2);
  pro_init(0) = w1;
  pro_init(1) = w2;
  arma::rowvec gamma_1(q);
  arma::rowvec gamma_2(q);
  arma::rowvec gamma_new(q);
  gamma_1 = Gamma.row(id1-1);
  gamma_2 = Gamma.row(id2-1);

  // Find those groups that are associated with the kth mixture
  int n1=getind(R,id1).n_cols;
  int n2=getind(R,id2).n_cols;
  int nid = n1 + n2;
  arma::rowvec ind1(m);
  arma::rowvec ind2(m);
  if(n1>0){
    ind1.set_size(n1);
    ind1 = getind(R,id1);
  }
  if(n2>0){
    ind2.set_size(n2);
    ind2 = getind(R,id2);
  }

  //re-arrange labels
  arma::rowvec R_new(m);
  R_new = R;
  for(int i=0; i<m; i++){
    if((R(i)>id1) & (R(i)<id2)){
      R_new(i)= R(i) -1;
    }
    if((R(i)>id2)){
      R_new(i)= R(i) -2;
    }
    if((R(i)==id1) || (R(i)==id2)){
      R_new(i) = Knew;
    }
  }

  //Compute the merge transformation using moment-matching principle
  double w1_new = w1 + w2;
  arma::rowvec w_new(Knew);
  w_new(Knew - 1) = w1_new;
  for(int i=0; i<Knew-1; i++){
    if(i<id1-1){
      w_new(i) = w(i);
    }
    if((i>=id1-1) & (i<id2-1)){
      w_new(i) = w(i+1);
    }
    if(i>=id2-1){
      w_new(i) = w(i+2);
    }
  }

  double alpha = w1/w1_new;
  arma::mat Gamma_new(Knew,q);
  gamma_new = (w1 * gamma_1 + w2 * gamma_2)/w1_new;
  Gamma_new = rbind0(removeiii(removeiii(Gamma,id1), id2-1),gamma_new);

  arma::rowvec beta(q);
  double ratio_propose = 0;
  for(int i=0; i<q; i++){
    beta(i) = (gamma_new(i) - gamma_1(i)) * std::sqrt(alpha/(1-alpha));
    if((beta(i)<1) & (beta(i)>0)){
      ratio_propose += R::dbeta(beta(i),2,2,true);
    }
  }
  //compute the proposal ratio
  arma::mat temp(1,1);
  if(n1>0){
    for(int i=0; i<n1; i++){
      pro = pro_init;
      ratio_propose += log(pro(0));
    }
  }
  if(n2>0){
    for(int i=0; i<n2; i++){
      pro = pro_init;
      ratio_propose += log(pro(1));
    }
  }

  //compute the log Jacobian of moment-matching transformation
  double log_Jacobian = log(w1_new) - (q/2.0) * log(alpha*(1-alpha));

  //compute the prior ratio
  double ratio_prior_w = (hyper_delta - 1 + n1) * log(w1) + (hyper_delta - 1 + n2) * log(w2) -
    (hyper_delta - 1 + nid) * log(w1_new) - log(R::beta(hyper_delta, Knew * hyper_delta));
  arma::mat C_new(Knew, Knew);
  C_new = updateC(Gamma_new,theta,tau);
  double ratio_prior_gamma = log(det_via_chol(C)) - log(det_via_chol(C_new));
  double ratio_prior = ratio_prior_w + ratio_prior_gamma;
  double death_ratio = log(K)*2 + ratio_propose  - ratio_prior - log_Jacobian;
  if(Knew ==1){
    death_ratio = death_ratio + log(2);
  }

  if(log(R::runif(0,1)) < death_ratio){
    K=Knew;
    Gamma = Gamma_new;
    R = R_new;
    C = C_new;
    w = w_new;
  }

  ans["K"] = K;
  ans["w"] = w;
  ans["Gamma"] = Gamma;
  ans["R"]= R;
  ans["C"] = C;
  return (ans);
}


List RJi_empty(arma::rowvec w,
               int K,
               arma::mat Gamma,

               arma::rowvec R,

               arma::mat C,
               double theta, double tau,
               int m,  int q,  double hyper_delta = 1){

  double merge_or_split = R::runif(0,1);
  if (((K == 1) || (merge_or_split < 0.5)) & (K<10)){
    return (
        Split_empty(w, K, Gamma, R,  C, theta, tau, m,  q,  hyper_delta)
    );
  }

  else{
    return (
        Merge_empty(w, K, Gamma, R, C, theta, tau, m,  q,  hyper_delta)
    );
  }
}
