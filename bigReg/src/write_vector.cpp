/*
  Date: 2016-07-23
  Author: Chibisi Chima-Okereke
*/
/* [[Rcpp::interfaces(r, cpp)]]*/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <fstream>

using namespace std;
using namespace arma;
using namespace Rcpp;

//' @title writes numeric vector to file
//' @description writes numeric vector to file
//' @param v numeric vector
//' @param filePath dependent variable
//' @export
// [[Rcpp::export]]
SEXP writeNumericVector(NumericVector v, std::string filePath){
	ofstream file;
	file.open(filePath.c_str(), ios::out | ios_base::binary);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.write(pointer, bytes);
	file.close();
	return wrap(0);
}

//' @title reads numeric vector to file
//' @description reads numeric vector to file
//' @param size the length of the numeric vector
//' @param filePath dependent variable
//' @export
// [[Rcpp::export]]
SEXP readNumericVector(int size, std::string filePath){
	ifstream file;
	file.open(filePath.c_str(), ios::in | ios::binary); 
	NumericVector v(size);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.read(pointer, bytes);
	file.close();
	return wrap(v);
}

//' @title creates factor from numeric vector and character vector as levels
//' @description The \code{CreateFactor} function creates a factor from a numeric 
//' 				vector and a character vector for levels
//' @param x numeric vector containing the numeric indices of the levels
//' @param levels character vector levels
//' @export
// [[Rcpp::export]]
SEXP CreateFactor(NumericVector x, CharacterVector levels)
{
  IntegerVector y = (IntegerVector)x;
  y.attr("levels") = levels;
  y.attr("class") = "factor";
  return wrap(y);
}

//' @title mySeq function to sequence integers
//' @description a function to create a sequence of integers
//' @param start integer from where sequence should start
//' @param end integer where sequence should end
//' @export
// [[Rcpp::export]]
IntegerVector mySeq(int start, int end){
	IntegerVector output(end - start + 1);
	for(int i = 0; i < output.size(); ++i){
		output[i] = start+i;
	}
	return output;
}

// TODO: Update with something better!
//' @title finds whether x is in y
//' @description finds whether x is in y
//' @param x item to be sought
//' @param y vector to be matched against
//' @export
// [[Rcpp::export]]
bool myIn(int x, NumericVector y){
	bool in = false;
	int i = 0;
	while(!in && i < y.size()){
		if(x == y[i]){
			in = true;
			break;
		}
		++i;
	}
	return in;
}

//' @title converts numeric vector to integer
//' @description converts numeric vector to integer
//' @param x numeric vector
//' @export
// [[Rcpp::export]]
IntegerVector asInteger(NumericVector x){
	IntegerVector y(x.size());
	for(int i = 0; i < x.size(); ++i){
		y[i] = x[i];
	}
	return y;
}

//' @title read data frame block from file
//' @description read data frame block from file
//' @param size number of elements in the block
//' @param filePath path to where the block is stored
//' @param df an empty list having the same number of elements as columns in the table
//' @param ncol number of columns in the dataframe block
//' @param factors list containing factors
//' @param factor_indices numeric vector containing the indicies that denote the factors
//' @export
// [[Rcpp::export]]
SEXP read_df_block(double size, std::string filePath, List df, 
				double ncol, List factors, NumericVector factor_indices){
	ifstream file;
	file.open(filePath.c_str(), ios::in | ios::binary); 
	arma::vec v(size);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.read(pointer, bytes);
	file.close();
	arma::vec y;
	NumericVector z;
	IntegerVector a;
	double nrow = size/ncol;
	int factor_index = 0;
	for(int i = 0; i < df.size(); ++i){
		y = v.subvec(i*nrow, (i + 1)*nrow - 1);
		z = wrap(y);
		if(myIn(i, factor_indices)){
			a = asInteger(z);
			a.attr("levels") = factors[factor_index];
			a.attr("class") = "factor";
			++factor_index;
			df[i] = a;
		}else{
			df[i] = z;
		}
	}
	df.attr("row.names") = mySeq(1, nrow);
	df.attr("class") = "data.frame";
	return wrap(df);
}

//' @title read multiple blocks of data frames from file
//' @description read multiple blocks of data frames from file
//' @param size number of elements in each block
//' @param filePaths path to where the blocks are stored
//' @param df an empty list having the same number of elements as columns in the table
//' @param ncols number of columns in the dataframe block
//' @param factors list containing factors
//' @param factor_indices numeric vector containing the indicies that denote the factors
//' @export
// [[Rcpp::export]]
SEXP read_df_blocks(NumericVector size, CharacterVector filePaths, List df, 
				double ncols, List factors, NumericVector factor_indices){

	ifstream file;
	std::string filePath;
	filePath = filePaths[0];
	file.open(filePath.c_str(), ios::in | ios::binary);
	arma::vec v(size[0]);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.read(pointer, bytes);
	file.close();
	arma::mat A = reshape(v, size[0]/ncols, ncols);

	for(int i = 1; i < filePaths.size(); ++i){
		ifstream file;
		filePath = filePaths[i];
		file.open(filePath.c_str(), ios::in | ios::binary);
		arma::vec v(size[i]);
		char* pointer = reinterpret_cast<char*>(&v[0]);
		size_t bytes = v.size() * sizeof(v[0]);
		file.read(pointer, bytes);
		file.close();
		// SLOW!!!
		A = join_cols(A, reshape(v, size[i]/ncols, ncols));
	}

	arma::vec y;
	NumericVector z;
	IntegerVector a;
	double nrow = A.n_rows;
	int factor_index = 0;
	for(int i = 0; i < df.size(); ++i){
		y = A.col(i);
		z = wrap(y);
		if(myIn(i, factor_indices)){
			a = asInteger(z);
			a.attr("levels") = factors[factor_index];
			a.attr("class") = "factor";
			++factor_index;
			df[i] = a;
		}else{
			df[i] = z;
		}
	}
	df.attr("row.names") = mySeq(1, nrow);
	df.attr("class") = "data.frame";
	return wrap(df);
}

//' @title writes numeric vector to file
//' @description writes numeric vector to file
//' @param v numeric vector to be written to file
//' @param filePath path to file where the numeric vector should be written
//' @export
// [[Rcpp::export]]
SEXP write_numeric_vector(arma::vec v, std::string filePath){ // "numericVector.bin"
	ofstream file;
	file.open(filePath.c_str(), ios::out | ios_base::binary);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.write(pointer, bytes);
	file.close();
	return wrap(0);
}


//' @title row binding for benchmarking ...
//' @description row binding for benchmarking
//' @param x first matrix to be bound together
//' @param y second matrix to be bound together
//' @export
// [[Rcpp::export]]
SEXP r_bind(arma::mat x, arma::mat y){
	return wrap(join_cols(x, y));
}

//' @title read matrix block from file
//' @description read matrix block from file
//' @param filePath path to file where matrix should be read from
//' @param size total number of elements to be read
//' @param ncol number of columns in the matrix
//' @export
// [[Rcpp::export]]
SEXP read_matrix_block(std::string filePath, double size, double ncol){
	ifstream file;
	file.open(filePath.c_str(), ios::in | ios::binary); 
	arma::vec v(size);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.read(pointer, bytes);
	file.close();
	NumericVector z = wrap(v);
	NumericVector dims(2);
	dims[0] = size/ncol;
	dims[1] = ncol;
	z.attr("class") = "matrix";
	z.attr("dim") = dims;
	return wrap(z);
}

//' @title read matrix blocks from file
//' @description read matrix blocks from file
//' @param filePaths file paths from where the matrix blocks will be read
//' @param size numeric vector containing the number of elements in each block
//' @param ncols number of columns in the matrix
//' @export
// [[Rcpp::export]]
SEXP read_matrix_blocks(CharacterVector filePaths, NumericVector size, double ncols){

	ifstream file;
	std::string filePath;
	filePath = filePaths[0];
	file.open(filePath.c_str(), ios::in | ios::binary);
	arma::vec v(size[0]);
	char* pointer = reinterpret_cast<char*>(&v[0]);
	size_t bytes = v.size() * sizeof(v[0]);
	file.read(pointer, bytes);
	file.close();
	arma::mat A = reshape(v, size[0]/ncols, ncols);

	for(int i = 1; i < filePaths.size(); ++i){
		ifstream file;
		filePath = filePaths[i];
		file.open(filePath.c_str(), ios::in | ios::binary);
		arma::vec v(size[i]);
		char* pointer = reinterpret_cast<char*>(&v[0]);
		size_t bytes = v.size() * sizeof(v[0]);
		file.read(pointer, bytes);
		file.close();
		// SLOW!!!
		A = join_cols(A, reshape(v, size[i]/ncols, ncols));
	}

	NumericMatrix B = wrap(A);
	NumericVector dims(2);
	dims[0] = A.n_rows;
	dims[1] = A.n_cols;
	B.attr("class") = "matrix";
	B.attr("dim") = dims;
	return wrap(B);
}
