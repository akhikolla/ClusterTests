#include <vector>

using namespace std;

double convertToDouble(const string& s);
vector< vector<double> > readCSV_to_vec(ifstream& file);
vector<double> colMeans(vector< vector<double> > x);
vector<double> colSds(vector< vector<double> > x);
double dist(vector<double> x, vector<double> y);
vector<double> quantile(vector<double> x, vector<double> prob);
double dinvgamma(double x, double shape, double scale, bool logp);
double dnormC(double x, double mu, double sig, bool logp);
double log_prior(vector<double> theta, int n_par);

