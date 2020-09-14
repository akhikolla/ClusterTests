#include <Rcpp.h>
using namespace Rcpp;

#include<cmath>

//' @title For all points in matrix 1, return the distance to and index of the
//' nearest point in matrix 2.
//'
//' @description
//' Find the shortest distance between each point in one data set and the points
//' in a second set.
//'
//' This function determines the distance between every point in data set 1 and
//' the points in data set 2. Unlike this function's naive counterpart,
//' find_min_dists, this function divides the PALM/iPALM localization data
//' into blocks, operates on the data in each block, and then performs linking
//' operations on neighboring blocks.
//'
//' @param mOne A numeric matrix where each row is a localization and each
//' column is a spatial axis.
//' @param mTwo A numeric matrix with the same number of columns as mOne.
//'
//' @author Zach Colburn
//'
//' @examples
//' \dontrun{
//' set.seed(10)
//'
//' mOne <- as.matrix(data.frame(
//' x = rnorm(10),
//' y = rbinom(10, 100, 0.5),
//' z = runif(10)
//' ))
//'
//' mTwo <- as.matrix(data.frame(
//' x = rnorm(20),
//' y = rbinom(20, 100, 0.5),
//' z = runif(20)
//' ))
//'
//' .find_min_dists_cpp(mOne, mTwo)
//'}
//'
//' @import Rcpp
//'
//' @useDynLib Bioi, .registration = TRUE
//[[Rcpp::export(.find_min_dists_cpp)]]
Rcpp::DataFrame find_min_dists_cpp(NumericMatrix mOne, NumericMatrix mTwo) {
  // Get number of points in each matrix
  //
  // mon = matrix one number of points
  // mtn = matrix two number of points
  int mon = mOne.nrow();
  int mtn = mTwo.nrow();

  // Get data dimensionality
  int dim = mOne.ncol();

  // Create matched index vector and min dist vector
  NumericVector matchedIndices(mon, 0.0);
  NumericVector minDists(mon, 0.0);



  // Iterate through data and update matched indices and min distances
  //
  // hp = home point (i.e. from matrix one)
  // fp = foreign point (i.e. from matrix two)
  for(int hp=0;hp<mon;hp++){
    // bi = best index (i.e. the index that minimizes the euclidean distance)
    // bsd = best squared distance
    //
    // The bsd is set to the value for the first check. This is done outside of
    // the main loop in order to get a baseline value or minimum distance. If a
    // value exceeds this at any point during the distance calculation, then it
    // is not the minimum. The value is updated during each loop of the main
    // calculation.
    int bi = 0;
    double bsd = 0.0;
    for(int i=0;i<dim;i++){
      bsd += pow(mOne(hp,i)-mTwo(0,i),2);
    }


    for(int fp=1;fp<mtn;fp++){
      // csd = current squared distance
      double csd = 0.0;
      for(int i=0;i<dim;i++){
        csd += pow(mOne(hp,i)-mTwo(fp,i),2);
        if(csd > bsd){
          break;
        }
      }
      if(csd<bsd){
        bsd = csd;
        bi = fp;
      }
    }
    matchedIndices(hp) = bi;
    minDists(hp) = sqrt(bsd);
  }

  for(int i=0;i<mon;i++){
    matchedIndices(i) += 1;
  }

  Rcpp::DataFrame output = Rcpp::DataFrame::create(
    Rcpp::Named("dist")=minDists,
    Rcpp::Named("index")=matchedIndices
  );

  return output;
}
