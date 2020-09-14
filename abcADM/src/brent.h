#include <vector>
namespace brent {

class func_base{
public:
  virtual double operator() (double, double[]) = 0;
};

double r8_abs ( double x );
double r8_epsilon ( );
double r8_max ( double x, double y );
double r8_sign ( double x );
void timestamp ( );
double zero ( double a, double b, double t, func_base& f );
void zero_rc ( double a, double b, double t, double &arg, int &status,
  double value );

// === simple wrapper functions
// === for convenience and/or compatibility
double glomin ( double a, double b, double c, double m, double e, double t,
  double f ( double x ), double &x );
double local_min ( double a, double b, double t, double f ( double x ),
  double &x );
double zero ( double a, double b, double t, double f ( double x, double params[] ), double params[] );

}
