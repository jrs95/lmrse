//Includes/namespaces
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

//' @title
//' robustseEigen
//' @description
//' robustseEigen fits cluster robust standard errors with C++Eigen across markers.
//' @param es a matrix of residuals.
//' @param xs a matrix of covariates.
//' @param cs a vector of clusters.
//' @param xxs a matrix of XtX.
//' @param covs an empty variance-covariance matrix.
//' @return A matrix of robust standard errors where the rows are the CpGs and the columns are the covariates including the intercept.
//' @author James R Staley <js16174@bristol.ac.uk>
//' @export
// [[Rcpp::export]]
NumericMatrix robustseEigen(NumericMatrix es, NumericMatrix xs, NumericMatrix cs, NumericMatrix xxs, NumericMatrix covs){
  const Eigen::Map<MatrixXd> e(as<Map<MatrixXd> >(es));
  const Eigen::Map<MatrixXd> x(as<Map<MatrixXd> >(xs));
  const Eigen::Map<MatrixXd> c(as<Map<MatrixXd> >(cs));
  const Eigen::Map<MatrixXd> xx(as<Map<MatrixXd> >(xxs));
  Map<MatrixXd> cov(as<Map<MatrixXd> >(covs));

  for(int j=0; j<e.cols(); ++j){

  int m = x.cols();
  int n = x.rows();
  MatrixXd Y(n,m) ;

  for(int i = 0; i < m; ++i) {
    Y.col(i) = e.col(j);
  }

  MatrixXd uj = c.adjoint() * x.cwiseProduct(Y);
  MatrixXd uu = (uj).adjoint() * uj;
  MatrixXd var = (xx * uu) * xx;
  for(int i = 0; i < m; ++i) {
    cov(j,i) = var(i,i);
  }
  }
  return wrap(cov);
}
