// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

double quad_tridiag(const arma::vec& x, const arma::vec& diag, double off_diag){
  double out = diag(0) * x(0) * x(0);
  
  for(arma::uword i = 1; i < x.n_elem; i++){
    out += diag(i) * x(i) * x(i) + 2.0 * off_diag * x(i) * x(i - 1);
  }
  
  return out;
}

double ldmvnorm_diag(const arma::vec& x, const arma::vec& mean, const arma::vec& diag){
  int n = x.n_elem;  // Dimensionality of the multivariate normal
  double log_det = arma::sum(arma::trunc_log(diag));  // Log of the determinant of diagonal covariance
  double log_density = - 0.5 * n * std::log(2 * M_PI) - 0.5 * log_det;
  
  // Compute the exponent term (squared Mahalanobis distance)
  arma::vec diff = x - mean;
  double exponent = - 0.5 * arma::sum(arma::square(diff) / diag);
  
  // Combine the terms to get the log density
  log_density += exponent;
  
  // Return the density in the original scale by exponentiating
  return log_density;
}

// Compute log-density for a multivariate normal with a tridiagonal precision matrix
// where all off-diagonal elements are equal

double ldmvnorm_tridiag_equal_offdiag(const arma::vec& x, const arma::vec& mean, const arma::vec& diag, double off_diag) {
  const arma::uword n = x.n_elem;
  
  // Check that diag has the correct dimensions
  if (diag.n_elem != n) {
    Rcpp::stop("Dimensions of 'diag' do not match with 'x' and 'mean'.");
  }
  
  // Calculate (x - mean)
  arma::vec diff = x - mean;
  
  // Compute log-determinant of tridiagonal matrix Q 
  
  arma::vec log_deltas(n, arma::fill::zeros);
  
  log_deltas(0) = arma::trunc_log(diag(0));
  
  for(arma::uword i = 1; i < n; i++) {
    log_deltas(i) = arma::trunc_log(diag(i)) + 
      arma::trunc_log(1.0 - (off_diag*off_diag)/(diag(i)*arma::trunc_exp(log_deltas(i - 1))));
  }
  
  double log_det = arma::sum(log_deltas);
  
  double exponent = diag(0) * diff(0) * diff(0);
  
  for(arma::uword i = 1; i < n; i++){
    exponent += diag(i) * diff(i) * diff(i) + 2.0 * off_diag * diff(i) * diff(i - 1); 
  }
  
  // Multivariate normal log-density
  double log_density = - 0.5 * (n * std::log(2 * M_PI) - log_det + exponent);
  
  return log_density;
}

List scaleandperiods(const arma::mat& data, bool scale) {
  // Get the dimensions of the data
  //int n = data.n_rows;
  int m = data.n_cols;
  
  // Initialize maxima and minima matrices
  arma::vec maxima(m);
  arma::vec minima(m);
  
  // Calculate maxima and minima for each column (time series)
  for (int i = 0; i < m; i++) {
    maxima(i) = data.col(i).max();
    minima(i) = data.col(i).min();
  }
  
  // Identify columns with constant time series
  arma::uvec cts = arma::find(maxima == minima);
  
  // If there are constant time series, remove them
  arma::mat mydata = data;
  if (!cts.is_empty()) {
    Rcpp::Rcout << "Removing series ";
    for (unsigned int i = 0; i < cts.n_elem; i++) {
      Rcpp::Rcout << (cts(i) + 1) << " ";
    }
    Rcpp::Rcout << "since they are constant." << std::endl;
    
    // Remove constant columns
    mydata.shed_cols(cts);
    m = mydata.n_cols;
    
    // Recompute maxima and minima for the remaining data
    maxima.resize(m);
    minima.resize(m);
    for (int i = 0; i < m; i++) {
      maxima(i) = mydata.col(i).max();
      minima(i) = mydata.col(i).min();
    }
  }
  
  // Scale data if required
  if (scale) {
    for (int j = 0; j < m; j++) {
      double range = maxima(j) - minima(j);
      if (range != 0) {
        mydata.col(j) = 1 + (1.0 / range) * (mydata.col(j) - maxima(j));
      }
    }
  }
  
  // Output the list
  return List::create(
    Named("mydata") = mydata,
    Named("cts") = cts
  );
}

arma::mat rmvnorm_arma_inv (int n, const arma::mat& precision, const arma::vec& location){
  
  int T = precision.n_rows;
  
  // arma::mat epsilon(T, n);
  // for (int i = 0; i < n; i++){
  //   epsilon.col(i) = as<arma::vec>(Rcpp::rnorm(T));
  // }
  
  arma::mat epsilon = arma::randn(T, n);
  
  arma::mat location_matrix(T, n, arma::fill::zeros);
  location_matrix.each_col() += location;
  
  arma::mat precision_chol_inv = trans(arma::inv(arma::trimatu(arma::chol(precision))));
  arma::mat draw = trans(precision_chol_inv) * (precision_chol_inv * location_matrix + epsilon);
  
  return draw;
}

arma::mat rmvnorm_arma_solve(int n, const arma::mat& precision, const arma::vec& location){
  
  int T = precision.n_rows;
  
  arma::mat epsilon = arma::randn(T, n);
  
  arma::mat location_matrix(T, n);
  location_matrix.each_col() = location;
  
  arma::mat precision_chol = arma::trimatu(arma::chol(precision));
  arma::mat draw = arma::solve(precision_chol, arma::solve(trans(precision_chol), location_matrix) + epsilon);
  
  return draw;
}

// super efficient draw from multivariate normal

List cholesky_tridiagonal(const arma::vec& omega_diag, const double& omega_offdiag) {
  const int T = omega_diag.n_elem - 1;
  arma::vec chol_diag(T + 1);
  arma::vec chol_offdiag(T + 1);
  chol_diag(0) = std::sqrt(omega_diag(0));
  for (int j = 1; j < T + 1; j++) {
    chol_offdiag(j - 1) = omega_offdiag / chol_diag(j-1);
    chol_diag(j) = std::sqrt(omega_diag(j) - chol_offdiag(j-1) * chol_offdiag(j-1));
  }
  return List::create(Named("chol_diag")=chol_diag, 
                      Named("chol_offdiag")=chol_offdiag);
}

arma::vec forward_algorithm(const arma::vec& chol_diag, const arma::vec& chol_offdiag, const arma::vec& covector) {
  const int T = chol_diag.n_elem - 1;
  arma:: vec htmp(T + 1);
  htmp(0) = covector(0) / chol_diag(0);
  for (int j = 1; j < T + 1; j++) {
    htmp(j) = (covector(j) - chol_offdiag(j - 1) * htmp(j - 1)) / chol_diag(j);
  }
  return htmp;
}

arma::vec backward_algorithm(const arma::vec& chol_diag, const arma::vec& chol_offdiag, const arma::vec& htmp) {
  const int T = chol_diag.size() - 1;
  arma::vec h(T + 1);
  h(T) = htmp(T) / chol_diag(T);
  for (int j = T - 1; j >= 0; j--) {
    h(j) = (htmp(j) - chol_offdiag(j) * h(j + 1)) / chol_diag(j);
  }
  return h;
}

arma::mat rmvnorm_arma_stochvol(int n, const arma::mat& precision, const arma::vec& location){
  
  int T = precision.n_rows;
  arma::vec precision_diag = precision.diag();
  double precision_offdiag = precision(1,0);
  
  List precision_chol = cholesky_tridiagonal(precision_diag, precision_offdiag);
  arma::vec aa = forward_algorithm(precision_chol["chol_diag"], precision_chol["chol_offdiag"], location);
  arma::mat draw(T, n);
  
  // arma::vec epsilon;
  arma::mat epsilon = arma::randn(T, n);
  for (int i = 0; i < n; i++){
    draw.col(i) = backward_algorithm(precision_chol["chol_diag"], precision_chol["chol_offdiag"], aa + epsilon.col(i));
  }
  
  return draw;
}

// [[Rcpp::export]]
List tseries_cpp(arma::mat data, arma::mat Z,
                 int runs, int burn, int thin,
                 arma::mat init_beta, 
                 arma::vec init_sig2beta, arma::vec init_sigma2, 
                 arma::vec init_phi, arma::vec init_tau2,
                 bool scale = false, int frequency = 365, int seasonfreq = 4, int seasondelay = 79, int deg = 2,
                 double a_eps = 1.0, double b_eps = 1.0, double a_beta = 1.0, double b_beta = 1.0, 
                 double a_tau = 1.0, double b_tau = 1.0, double a_phi = 1.0, double b_phi = 1.0, 
                 double eta_phi = 1.0, 
                 bool PPM = false) {
  
  //// CONSTRUCTION OF THE DESIGN MATRICES ////
  
  List data_scaled = scaleandperiods(data, scale);
  
  // Extract elements from the returned list
  arma::mat y = data_scaled["mydata"];
  // arma::mat y = data;
  arma::uvec cts = data_scaled["cts"];
  
  // Remove rows from locations based on cts
  if (!cts.is_empty()) {
    // locations.shed_rows(cts);
    init_beta.shed_cols(cts);
    init_phi.shed_rows(cts);
    init_tau2.shed_rows(cts);
    init_sigma2.shed_rows(cts);
  }
  
  // Number of periods in the time series
  int T = y.n_rows;
  // Number of time series in the data
  int n = y.n_cols;
  // Number of columns of Z
  int p = Z.n_cols;

  // Initializing Sigma_n matrices
  arma::vec sigma2 = init_sigma2; 
  arma::cube Sigma_n(T, T, n, arma::fill::zeros);
  arma::cube invSigma_n(T, T, n, arma::fill::zeros);
  
  for (int i = 0; i < n; i++) {
    Sigma_n.slice(i) = sigma2(i) * arma::eye(T, T);
    invSigma_n.slice(i) = 1.0 / sigma2(i) * arma::eye(T, T);
  }
  
  arma::vec phi = init_phi;      
  arma::vec tau2 = init_tau2;
  
  // Initialize gamma as a concatenation of rho and sig2the
  arma::mat gamma(2, n);
  gamma.row(0) = phi.t();
  gamma.row(1) = tau2.t();
  
  arma::cube Phiinv(T, T, n, arma::fill::zeros);
  arma::cube Rinv(T, T, n, arma::fill::zeros);
  
  for (int i = 0; i < n; i++) {
    Phiinv.slice(i).diag(0).fill(1.0 + phi(i) * phi(i));
    Phiinv.slice(i).diag(1).fill(-phi(i));
    Phiinv.slice(i).diag(-1).fill(-phi(i));
    Phiinv(0, 0, i) = 1.0;
    Phiinv(T - 1, T - 1, i) = 1.0;
    Rinv.slice(i) = 1.0 / tau2(i) * Phiinv.slice(i);
  }
  
  // Initializing sig2beta and variance-covariance matrices
  arma::vec sig2beta = init_sig2beta;
  arma::mat sigmabeta = arma::diagmat(sig2beta);
  arma::mat invsigmabeta = arma::diagmat(1.0 / sig2beta);
  
  arma::mat beta = init_beta; 
  arma::mat loc = y - Z * beta;
  
  arma::cube Q_n(T, T, n, arma::fill::zeros);
  arma::cube Vwinv_n(T, T, n, arma::fill::zeros);
  
  arma::mat w(T, n, arma::fill::zeros); 
  
  // Fill matrices and w for each i
  
  for (int i = 0; i < n; i++) {
    
    // Recursive generation of AR(1) samples
    w(0, i) = arma::randn<double>(arma::distr_param(0.0, sqrt(tau2(i))));  // First element
    for (int t = 1; t < T; t++) {
      w(t, i) = phi(i) * w(t - 1, i) + arma::randn<double>(arma::distr_param(0.0, sqrt(tau2(i))));
    }
    
    Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv.slice(i), arma::eye(T, T));
    Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv.slice(i); 
  }
  
  // Determine CL based on thinning
  int CL;
  if (thin == 0) {
    CL = std::floor(runs - burn);
  } else {
    CL = std::floor((runs - burn) / thin);
  }
  
  // Memory allocation for matrices and arrays
  arma::cube betasample(p, n, CL, arma::fill::zeros);   // beta samples
  arma::cube wsample(T, n, CL, arma::fill::zeros);      // w samples
  arma::mat sigma2sample(CL, n, arma::fill::zeros);     // sigma2 posterior samples
  arma::mat sig2betasample(CL, p, arma::fill::zeros);   // sig2beta posterior samples
  arma::mat phisample(CL, n, arma::fill::zeros);        // phi posterior samples
  arma::mat tau2sample(CL, n, arma::fill::zeros);       // tau2 posterior samples
  arma::vec loglik(CL, arma::fill::zeros);
  
  int acc_phi = 0;
  
  for (int iter = 0; iter < runs; iter++) {
    
    double progress = static_cast<double>(iter * 100) / runs;
    if (progress >= 10.0 && progress <= 100.0 && fmod(progress, 10.0) == 0.0) {
      Rcout << "### Progress: " << progress << " %" << std::endl;
      Rcout << " " << std::endl;
    }
    
    /////////////////
    // beta update //
    /////////////////
    
    for (int i = 0; i < n; i++) {
      
      // Covariance matrix for beta update
      arma::mat Vbetainv = Z.t() * Q_n.slice(i) * Z + invsigmabeta;
      arma::mat Vbeta = arma::inv_sympd(Vbetainv);
      
      // Mean vector for beta update
      arma::vec mubeta = Vbeta * Z.t() * Q_n.slice(i) * y.col(i);
      
      // Sample new alpha from multivariate normal distribution
      beta.col(i) = arma::mvnrnd(mubeta, Vbeta, 1);
      
    }
    
    loc = y - Z * beta;
    
    //////////////
    // w update //
    //////////////
    
    for (int i = 0; i < n; i++){
      w.col(i) = rmvnorm_arma_stochvol(1, Vwinv_n.slice(i), 1 / sigma2(i) * loc.col(i));
    }
    
    ///////////////////
    // update sigma2 //
    ///////////////////
    
    arma::mat residuals = loc - w;
    
    // Compute sigma2 using inverse gamma
    for (int i = 0; i < n; i++) {
      double shape = a_eps + T / 2.0;
      double scale = b_eps + arma::sum(arma::square(residuals.col(i))) / 2.0;
      
      sigma2(i) = 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / scale));
      
      Sigma_n.slice(i) = sigma2(i) * arma::eye(T, T);
      invSigma_n.slice(i) = 1.0 / sigma2(i) * arma::eye(T, T);
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv.slice(i), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv.slice(i);
    }
    
    ///////////////////////
    // update sigma2beta //
    ///////////////////////
    
    // Compute `sig2beta` using inverse gamma
    for (int j = 0; j < p; j++) {
      double shape = a_beta + n / 2.0;
      double scale = b_beta + arma::sum(arma::square(beta.row(j))) / 2.0;
      
      sig2beta(j) = 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / scale));
      
    }
    
    // Construct sigmaalpha and invsigmaalpha as diagonal matrices
    arma::mat invsigmabeta = arma::diagmat(1.0 / sig2beta);
    
    ////////////////////
    // update phistar //
    ////////////////////
    
    arma::vec trans_phi = arma::trunc_log(phi) - arma::log1p(-phi);
    
    arma::vec trans_phi_prop = arma::mvnrnd(trans_phi, eta_phi * arma::eye(n, n), 1);
    
    arma::vec phi_prop = arma::trunc_exp(trans_phi_prop) / (1 + arma::trunc_exp(trans_phi_prop));
    
    arma::cube Phiinv_prop(T, T, n, arma::fill::zeros);
    arma::cube Rinv_prop(T, T, n, arma::fill::zeros);
    
    arma::vec log_target = (a_phi - 1) * arma::trunc_log(phi) +
      (b_phi - 1) * arma::log1p(-phi) +
      (trans_phi - 2 * arma::log1p(arma::trunc_exp(trans_phi)));
    
    arma::vec log_target_prop = (a_phi - 1) * arma::trunc_log(phi_prop) +
      (b_phi - 1) * arma::log1p(-phi_prop) +
      (trans_phi_prop - 2 * arma::log1p(arma::trunc_exp(trans_phi_prop)));
    
    for(int i = 0; i < n; i++){
      
      Phiinv_prop.slice(i).diag(0).fill(1 + phi_prop(i) * phi_prop(i));
      Phiinv_prop.slice(i).diag(1).fill(-phi_prop(i));
      Phiinv_prop.slice(i).diag(-1).fill(-phi_prop(i));
      Phiinv_prop(0, 0, i) = 1.0;
      Phiinv_prop(T - 1, T - 1, i) = 1.0;
      
      Rinv_prop.slice(i) = 1.0 / tau2(i) * Phiinv_prop.slice(i);
      
      log_target_prop(i) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_prop.slice(i).diag(0), Rinv_prop(0, 1, i));
      log_target(i) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv.slice(i).diag(0), Rinv(0, 1, i));
      
      double log_ratio = log_target_prop(i) - log_target(i);
      
      double log_alpha = std::min(0.0, log_ratio);
      
      double u = arma::randu();
      
      if(log(u) < log_alpha){
        
        acc_phi += 1;
        phi(i) = phi_prop(i);
        Phiinv.slice(i) = Phiinv_prop.slice(i);
        Rinv.slice(i) = Rinv_prop.slice(i);
        
      }
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv.slice(i), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv.slice(i);
      
    }
    
    /////////////////////
    // update tau2star //
    /////////////////////
    
    for(int i = 0; i < n; i++){
      double shape = a_tau + T / 2.0;
      double scale = b_tau + quad_tridiag(w.col(i), Phiinv.slice(i).diag(0), Phiinv(0, 1, i)) / 2.0;
      
      tau2(i) = 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / scale));
      
      Rinv.slice(i) = 1.0 / tau2(i) * Phiinv.slice(i);
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv.slice(i), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv.slice(i);
    }
    
    /////////////////////
    // saving iterates //
    /////////////////////
    
    if((fmod(iter, thin) == 0.0) & (iter >= burn)){
      betasample.slice(iter - burn) = beta;
      wsample.slice(iter - burn) = w;
      sigma2sample.row(iter - burn) = sigma2.t();
      sig2betasample.row(iter - burn) = sig2beta.t();
      phisample.row(iter - burn) = phi.t();
      tau2sample.row(iter - burn) = tau2.t();
      for (int i = 0; i < n; i++) {
        loglik(iter - burn) += ldmvnorm_diag(y.col(i), Z*beta.col(i) + w.col(i), sigma2(i) * arma::ones(T));
      }
    }
    
  }
  
  //return output;
  return Rcpp::List::create(
    Rcpp::Named("beta") = betasample,
    Rcpp::Named("w") = wsample,
    Rcpp::Named("sigma2") = sigma2sample,
    Rcpp::Named("sigma2beta") = sig2betasample,
    Rcpp::Named("phi") = phisample,
    Rcpp::Named("tau2") = tau2sample,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("rem_obs") = cts + 1
  );
  
}