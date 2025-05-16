// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

int sample_cpp(arma::uword n, int size, const arma::vec& probs, bool replace = true) {
  if (probs.n_elem != n) {
    Rcpp::stop("Length of probabilities vector must match n.");
  }
  // Convert arma::vec to Rcpp::NumericVector
  Rcpp::NumericVector rcpp_probs(probs.begin(), probs.end());
  
  // Perform sampling
  return Rcpp::sample(n, size, replace, rcpp_probs)[0];
}

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
  // return List::create(_["chol_diag"]=chol_diag, _["chol_offdiag"]=chol_offdiag);
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

// Similarity 1

double Similarity1(const arma::vec& s1, const arma::vec& s2, double epsilon, int dim, int lg){
  
  double cent1 = 0.0; 
  double cent2 = 0.0;
  double out = 0.0;
  double sdist, maxdist, dist;
  
  for(int ii = 0; ii < dim; ii++){
    cent1 = cent1 + s1(ii)/(double) dim;
    cent2 = cent2 + s2(ii)/(double) dim;
  }
  
  maxdist = sqrt((s1(0) - cent1)*(s1(0) - cent1) +
    (s2(0) - cent2)*(s2(0) - cent2));
  
  sdist = 0;
  
  for(int ii = 0; ii < dim; ii++){
    dist = sqrt((s1(ii) - cent1)*(s1(ii) - cent1) +
      (s2(ii) - cent2)*(s2(ii) - cent2));
    sdist = sdist + dist;
    if(maxdist < dist) {
      maxdist = dist;
    } 
  }
  
  if (sdist >= 1) {
    out = -lgamma(sdist * epsilon);
  } else if (sdist != 0) {
    out = -arma::trunc_log(sdist);
  } else {
    out = log(1);
  }
  if (!lg) {
    out = arma::trunc_exp(out);
  } 
  return(out);
}

// Similarity 2

double Similarity2(const arma::vec& s1, const arma::vec& s2, double a, int dim, int lg){
  
  double dist = 0.0;
  double out = 0.0;
  
  for(int i = 0; i < dim; i++){
    for(int ii = 0; ii < dim; ++ii){
      dist = sqrt((s1(i) - s1(ii))*(s1(i) - s1(ii)) +
        (s2(i) - s2(ii))*(s2(i) - s2(ii)));
      if (dist < a) {
        out = 1.0;
      } else {
        out = 0.0;
        break;
      }
    }
    if (out < 1) {
      break;
    } 
  }
  
  if(lg) { 
    out = arma::trunc_log(out);
  } 
  
  return(out);
  
}


// Similarity 3 and 4

double G2a(double a, int lg){
  
  double out;
  
  out = log(M_PI) + lgamma(a) + lgamma(a - 0.5);
  
  if(!lg) {
    out = arma::trunc_exp(out);
  } 
  
  return(out);
}

double Similarity3_4(const arma::vec& s1, const arma::vec& s2, const arma::vec& mu0, double k0, double v0, const arma::mat& L0,
                   int dim, int Similarity, int lg){
  
  
  double kn, knn, vn, vnn, sbar1 = 0.0, sbar2 = 0.0, dL0, dLn, dLnn;
  double out = 0.0;
  double s_sbar1 = 0.0;
  double s_sbar2 = 0.0;
  double sbar_mu01, sbar_mu02, sbar_mun1, sbar_mun2, mun1, mun2;
  double Vs1 = 0.0;
  double Vs2 = 0.0;
  double Vs3 = 0.0;
  double Vs4 = 0.0;
  double Vsbarmu01, Vsbarmu02, Vsbarmu03, Vsbarmu04;
  double Vsbarmun1, Vsbarmun2, Vsbarmun3, Vsbarmun4;
  double Ln1, Ln2, Ln3, Ln4, Lnn1, Lnn2, Lnn3, Lnn4;
  
  for(int ii = 0; ii < dim; ii++){
    
    sbar1 = sbar1 + s1(ii)/(double) dim;
    sbar2 = sbar2 + s2(ii)/(double) dim;
    
  }
  
  for(int ii = 0; ii < dim; ii++){
    
    s_sbar1 = s1(ii) - sbar1;
    s_sbar2 = s2(ii) - sbar2;
    
    
    Vs1 = Vs1 + s_sbar1 * s_sbar1;
    Vs2 = Vs2 + s_sbar1 * s_sbar2;
    Vs3 = Vs3 + s_sbar2 * s_sbar1;
    Vs4 = Vs4 + s_sbar2 * s_sbar2;
    
  }
  
  kn = k0 + dim; vn = v0 + dim;
  knn = kn + dim; vnn = vn + dim;
  
  mun1 = k0 / (k0 + dim) * mu0(0) + dim / (k0 + dim) * sbar1;
  mun2 = k0 / (k0 + dim) * mu0(1) + dim / (k0 + dim) * sbar2;
  
  sbar_mu01 = sbar1 - mu0(0);
  sbar_mu02 = sbar2 - mu0(1);
  sbar_mun1 = sbar1 - mun1;
  sbar_mun2 = sbar2 - mun2;
  
  Vsbarmu01 = sbar_mu01 * sbar_mu01, Vsbarmun1 = sbar_mun1 * sbar_mun1;
  Vsbarmu02 = sbar_mu01 * sbar_mu02, Vsbarmun2 = sbar_mun1 * sbar_mun2;
  Vsbarmu03 = sbar_mu02 * sbar_mu01, Vsbarmun3 = sbar_mun2 * sbar_mun1;
  Vsbarmu04 = sbar_mu02 * sbar_mu02, Vsbarmun4 = sbar_mun2 * sbar_mun2;
  
  Ln1 = L0(0, 0) + Vs1 + k0 * dim / (k0 + dim) * Vsbarmu01;
  Ln2 = L0(0, 1) + Vs2 + k0 * dim / (k0 + dim) * Vsbarmu02;
  Ln3 = L0(1, 0) + Vs3 + k0 * dim / (k0 + dim) * Vsbarmu03;
  Ln4 = L0(1, 1) + Vs4 + k0 * dim / (k0 + dim) * Vsbarmu04;
  
  Lnn1 = Ln1 + Vs1 + kn * dim / (kn + dim) * Vsbarmun1;
  Lnn2 = Ln2 + Vs2 + kn * dim / (kn + dim) * Vsbarmun2;
  Lnn3 = Ln3 + Vs3 + kn * dim / (kn + dim) * Vsbarmun3;
  Lnn4 = Ln4 + Vs4 + kn * dim / (kn + dim) * Vsbarmun4;
  
  dL0 = L0(0, 0) * L0(1, 1) - L0(0, 1) * L0(1, 0);
  dLn = Ln1 * Ln4 - Ln2 * Ln3;
  dLnn = Lnn1 * Lnn4 - Lnn2 * Lnn3;
  
  if(Similarity == 3){
    out = - dim * log(M_PI) +
      (G2a(0.5 * vn, 1) - G2a(0.5 * v0, 1)) +
      (0.5 * v0 * arma::trunc_log(dL0) - 0.5 * vn * arma::trunc_log(dLn)) +
      (arma::trunc_log(k0) - arma::trunc_log(kn));
  }
  if(Similarity == 4){
    out = - dim * log(M_PI) +
      (G2a(0.5 * vnn, 1) - G2a(0.5 * vn, 1)) +
      (0.5 * vn * arma::trunc_log(dLn) - 0.5 * vnn * arma::trunc_log(dLnn)) +
      (arma::trunc_log(kn) - arma::trunc_log(knn));
  }
  if(!lg) {
    out = arma::trunc_exp(out);  
  }
  
  return(out);
  
}

List unique_columns_info(const arma::mat& gamma) {
  int n = gamma.n_cols; // Number of columns in gamma
  
  // Containers for unique column indices, cardinalities, and labels
  arma::uvec unique_indices;
  arma::uvec cardinalities; // Cardinality of each unique column
  arma::uvec labels(n, arma::fill::zeros); // Group labels for each column
  
  // Iterate over all columns to find unique ones
  for (int j = 0; j < n; j++) {
    // Extract the current column
    arma::vec current_col = gamma.col(j);
    bool is_unique = true;
    
    // Check if this column matches any of the already identified unique columns
    for (size_t k = 0; k < unique_indices.n_elem; k++) {
      if (arma::all(gamma.col(unique_indices(k)) == current_col)) {
        is_unique = false;
        labels(j) = k + 1; // Assign the group label (1-based index)
        ++cardinalities(k); // Increment the cardinality of the matched group
        break;
      }
    }
    
    // If the column is unique, add it to the list of unique columns
    if (is_unique) {
      unique_indices.insert_rows(unique_indices.n_elem, arma::uvec({(unsigned int)j}));
      cardinalities.insert_rows(cardinalities.n_elem, arma::uvec({1})); // Initialize cardinality to 1
      labels(j) = unique_indices.n_elem; // Assign a new group label
    }
  }
  
  // Extract the unique columns using the indices
  arma::mat unique_columns = gamma.cols(unique_indices);
  
  // Return the results as a list
  return List::create(
    Named("unique_columns") = unique_columns,
    Named("cardinalities") = cardinalities,
    Named("labels") = labels
  );
}

// [[Rcpp::export]]
List tseriesclust_similarity_cpp(arma::mat data, arma::mat locations, arma::mat Z,
                                 int runs, int burn, int thin,
                                 int similarity, 
                                 arma::mat init_beta, 
                                 arma::vec init_sig2beta, arma::vec init_sigma2, 
                                 arma::vec init_phi, arma::vec init_tau2, double init_M,
                                 double k0 = 1.0, double v0 = 2.0, 
                                 double eC1 = 1.0, double aC2 = 1.5, 
                                 bool scale = false, int frequency = 365, int seasonfreq = 4, int seasondelay = 79, int deg = 2,
                                 double a_eps = 1.0, double b_eps = 1.0, double a_beta = 1.0, double b_beta = 1.0, 
                                 double a_tau = 1.0, double b_tau = 1.0, double a_phi = 1.0, double b_phi = 1.0, 
                                 double eta_phi = 1.0, 
                                 bool PPM = false, bool clust = true, int CC = 3) {
  
  // Compute column means of the locations matrix
  arma::vec mu0 = trans(arma::mean(locations, 0)); // Column means
  
  // Create an identity matrix with dimensions based on the number of columns in locations
  arma::mat L0 = arma::eye(locations.n_cols, locations.n_cols);
  
  //// CONSTRUCTION OF THE DESIGN MATRICES ////
  
  List data_scaled = scaleandperiods(data, scale);
  
  // Extract elements from the returned list
  arma::mat y = data_scaled["mydata"];
  // arma::mat y = data;
  arma::uvec cts = data_scaled["cts"];
  
  // Remove rows from locations based on cts
  if (!cts.is_empty()) {
    locations.shed_rows(cts);
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
  
  // Extract latitude and longitude columns
  arma::vec s1 = locations.col(1);
  arma::vec s2 = locations.col(0);
  
  arma::vec s1o(n, arma::fill::zeros);
  arma::vec s2o(n, arma::fill::zeros);
  
  arma::vec s1n(n, arma::fill::zeros);
  arma::vec s2n(n, arma::fill::zeros);
  
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
  
  double Mdp = init_M;
  
  // Initialize gamma as a concatenation of rho and sig2the
  arma::mat gamma(2, n);
  gamma.row(0) = phi.t();
  gamma.row(1) = tau2.t();
  
  // clustering variables
  
  List summary_gamma = unique_columns_info(gamma);
  
  arma::mat gamma_star = summary_gamma["unique_columns"];
  arma::vec phi_star = gamma_star.row(0).t();
  arma::vec tau2_star = gamma_star.row(1).t();
  
  int k_current = gamma_star.n_cols; 
  arma::ivec nj_current = summary_gamma["cardinalities"];
  arma::ivec clu_current = summary_gamma["labels"];
  
  arma::cube Phiinv_star(T, T, k_current, arma::fill::zeros);
  arma::cube Rinv_star(T, T, k_current, arma::fill::zeros);
  
  for (int k = 0; k < k_current; k++) {
    Phiinv_star.slice(k).diag(0).fill(1.0 + phi_star(k) * phi_star(k));
    Phiinv_star.slice(k).diag(1).fill(-phi_star(k));
    Phiinv_star.slice(k).diag(-1).fill(-phi_star(k));
    Phiinv_star(0, 0, k) = 1.0;
    Phiinv_star(T - 1, T - 1, k) = 1.0;
    Rinv_star.slice(k) = 1.0 / tau2_star(k) * Phiinv_star.slice(k);
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
    
    Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv_star.slice(clu_current(i) - 1), arma::eye(T, T));
    Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv_star.slice(clu_current(i) - 1); 
  }
  
  
  double lCo = 0.0;
  double lCn = 0.0;
  double lCn_1 = 0.0;
  
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
  arma::imat clustmemory(CL, n, arma::fill::zeros);
  arma::ivec k_clust(CL, arma::fill::zeros);
  arma::vec M_dp(CL, arma::fill::zeros);
  arma::vec loglik(CL, arma::fill::zeros);
  
  // Declare the  output containers
  List nj(CL);
  List phistarsample(CL);
  List tau2starsample(CL);
  
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
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv_star.slice(clu_current(i) - 1), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv_star.slice(clu_current(i) - 1);
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
    
    //////////
    // sPPM //
    //////////
    
    if (clust == true) {
      
      for (int i = 0; i < n; i++) {
        
        arma::vec phi_empty = Rcpp::rbeta(CC, a_phi, b_phi);
        arma::vec tau2_empty = 1.0 / arma::randg(CC, arma::distr_param(a_tau, 1.0 / b_tau));
        
        arma::cube Rinv_empty(T, T, CC, arma::fill::zeros);
        arma::cube Phiinv_empty(T, T, CC, arma::fill::zeros);
        
        for (int cc = 0; cc < CC; cc++) {
          
          Phiinv_empty.slice(cc).diag(0).fill(1 + phi_empty(cc) * phi_empty(cc));
          Phiinv_empty.slice(cc).diag(1).fill(-phi_empty(cc));
          Phiinv_empty.slice(cc).diag(-1).fill(-phi_empty(cc));
          Phiinv_empty(0, 0, cc) = 1.0;
          Phiinv_empty(T - 1, T - 1, cc) = 1.0;
          
          Rinv_empty.slice(cc) = 1.0 / tau2_empty(cc) * Phiinv_empty.slice(cc);
        }
        
        int zi = clu_current(i);
        
        if (nj_current(zi - 1) > 1){
          
          nj_current(zi - 1) -= 1;
          
          arma::vec wei_filled(k_current, arma::fill::zeros);
          
          for (int k = 0; k < k_current; k++) {
            
            ////////////////
            // SIMILARITY //
            ////////////////
            
            int indx_k = 0;
            
            for(int ii = 0; ii < n; ii++){
              
              if((clu_current(ii) == k + 1) & (i != ii)){
                
                s1o(indx_k) = s1(ii);
                s2o(indx_k) = s2(ii);
                
                s1n(indx_k) = s1(ii);
                s2n(indx_k) = s2(ii);
                
                indx_k = indx_k + 1;
                
              }
              
              if(i == ii){
                
                s1n(nj_current(k)) = s1(ii);
                s2n(nj_current(k)) = s2(ii);
                
              }
              
            }
            
            if(similarity==1){
              lCo = Similarity1(s1o, s2o, eC1, nj_current(k), 1);
              lCn = Similarity1(s1n, s2n, eC1, nj_current(k) + 1, 1);
            }
            
            if(similarity==2){
              lCo = Similarity2(s1o, s2o, aC2, nj_current(k), 1);
              lCn = Similarity2(s1n, s2n, aC2, nj_current(k) + 1, 1);
            }
            
            if(similarity==3){
              lCo = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, nj_current(k), 3, 1);
              lCn = Similarity3_4(s1n, s2n, mu0, k0, v0, L0, nj_current(k) + 1, 3, 1);
            }
            
            if(similarity==4){
              lCo = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, nj_current(k), 4, 1);
              lCn = Similarity3_4(s1n, s2n, mu0, k0, v0, L0, nj_current(k) + 1, 4, 1);
            }
            
            //////////////////////////////////////////////////////////////////////
            
            wei_filled(k) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_star.slice(k).diag(0), Rinv_star(0, 1, k)) +
              lCn - lCo + // Spatial part of similarity function
              log(nj_current(k));
            
          }
          
          arma::vec wei_empty(CC, arma::fill::zeros);
          
          for (int cc = 0; cc < CC; cc++) {
            
            s1o(0) = s1(i);
            s2o(0) = s2(i);
            
            if(similarity==1) lCn_1 = log(1); // This is because set the similarity to M with singletons
            if(similarity==2) lCn_1 = log(1);
            if(similarity==3) lCn_1 = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, 1, 3, 1);
            if(similarity==4) lCn_1 = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, 1, 4, 1);
            
            wei_empty(cc) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_empty.slice(cc).diag(0), Rinv_empty(0, 1, cc)) +
              lCn_1 + // Spatial part of similarity function
              log(Mdp) - log(CC);
            
          }
          
          arma::vec wei = arma::join_cols(wei_filled, wei_empty);
          
          // Trick to avoid all zero in the weights
          double max_wei = arma::max(wei);
          wei = wei - max_wei;
          wei = arma::trunc_exp(wei)/arma::sum(arma::trunc_exp(wei));
          
          int newci = sample_cpp(k_current + CC, 1, wei, true);
          
          // Check if newci is less than or equal to k_current
          if (newci <= k_current) {
            // Point 2a, page 345 of Favaro and Teh, second column
            clu_current(i) = newci;          // Update cluster assignment
            nj_current(newci - 1) += 1;      // Increment nj_current for the selected cluster
          } else {
            int id_empty = newci - k_current;
            
            // Point 2b  page 345 Favaro and Teh, second column
            clu_current(i) = k_current + 1;
            k_current += 1;
            
            nj_current.resize(k_current); // Resize nj_current to accommodate the new cluster
            nj_current(k_current - 1) = 1;
            
            phi_star.resize(k_current);
            phi_star(k_current - 1) = phi_empty(id_empty - 1);
            
            tau2_star.resize(k_current);
            tau2_star(k_current - 1) = tau2_empty(id_empty - 1);
            
            Phiinv_star = arma::join_slices(Phiinv_star, Phiinv_empty.slice(id_empty - 1));
            Rinv_star = arma::join_slices(Rinv_star, Rinv_empty.slice(id_empty - 1));
            
            phi_empty(id_empty - 1) = Rcpp::rbeta(1, a_phi, b_phi)[0];
            tau2_empty(id_empty - 1) = 1.0 / arma::randg<double>(arma::distr_param(a_tau, 1.0 / b_tau));
            
            Phiinv_empty.slice(id_empty - 1).diag(0).fill(1 + phi_empty(id_empty - 1) * phi_empty(id_empty - 1));
            Phiinv_empty.slice(id_empty - 1).diag(1).fill(-phi_empty(id_empty - 1));
            Phiinv_empty.slice(id_empty - 1).diag(-1).fill(-phi_empty(id_empty - 1));
            Phiinv_empty(0, 0, id_empty - 1) = 1.0;
            Phiinv_empty(T - 1, T - 1, id_empty - 1) = 1.0;
            
            Rinv_empty.slice(id_empty - 1) = 1.0 / tau2_empty(id_empty - 1) * Phiinv_empty.slice(id_empty - 1);
            
          }
          
        } else {
          // We are in a case where nj=1
          
          // Operations described on the paragraph starting
          // with On the other hand of  page 345 Favaro and Teh, second column
          
          int id_empty = Rcpp::sample(CC, 1, true)[0];
          
          phi_empty(id_empty - 1) = phi_star(zi - 1);
          tau2_empty(id_empty - 1) = tau2_star(zi - 1);
          Phiinv_empty.slice(id_empty - 1) = Phiinv_star.slice(zi - 1);
          Rinv_empty.slice(id_empty - 1) = Rinv_star.slice(zi - 1);
          
          // this is complicated I have to adjust the current state
          nj_current.shed_row(zi - 1);
          phi_star.shed_row(zi - 1);
          tau2_star.shed_row(zi - 1);
          Phiinv_star.shed_slice(zi - 1);
          Rinv_star.shed_slice(zi - 1);
          
          k_current -= 1;
          for(int ii = 0; ii < n; ii++){
            if(clu_current(ii) > zi){ // the cluster with label greater than the one i'm discarding have to decrease by one
              clu_current(ii) -= 1;
            }
          }
          
          arma::vec wei_filled(k_current, arma::fill::zeros);
          
          for (int k = 0; k < k_current; k++) {
            
            int indx_k = 0;
            
            for(int ii = 0; ii < n; ii++){
              
              if((clu_current(ii) == k + 1) & (i != ii)){
                
                s1o(indx_k) = s1(ii);
                s2o(indx_k) = s2(ii);
                
                s1n(indx_k) = s1(ii);
                s2n(indx_k) = s2(ii);
                
                indx_k = indx_k + 1;
                
              }
              
              if(i == ii){
                
                s1n(nj_current(k)) = s1(ii);
                s2n(nj_current(k)) = s2(ii);
                
              }
              
            }
            
            if(similarity==1){
              lCo = Similarity1(s1o, s2o, eC1, nj_current(k), 1);
              lCn = Similarity1(s1n, s2n, eC1, nj_current(k) + 1, 1);
            }
            
            if(similarity==2){
              lCo = Similarity2(s1o, s2o, aC2, nj_current(k), 1);
              lCn = Similarity2(s1n, s2n, aC2, nj_current(k) + 1, 1);
            }
            
            if(similarity==3){
              lCo = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, nj_current(k), 3, 1);
              lCn = Similarity3_4(s1n, s2n, mu0, k0, v0, L0, nj_current(k) + 1, 3, 1);
            }
            
            if(similarity==4){
              lCo = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, nj_current(k), 4, 1);
              lCn = Similarity3_4(s1n, s2n, mu0, k0, v0, L0, nj_current(k) + 1, 4, 1);
            }
            
            wei_filled(k) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_star.slice(k).diag(0), Rinv_star(0, 1, k)) +
              lCn - lCo + // Spatial part of similarity function
              log(nj_current(k));
            
          }
          
          arma::vec wei_empty(CC, arma::fill::zeros);
          
          for (int cc = 0; cc < CC; cc++) {
            
            s1o(0) = s1(i);
            s2o(0) = s2(i);
            
            if(similarity==1) lCn_1 = log(1); // This is because set the similarity to M with singletons
            if(similarity==2) lCn_1 = log(1);
            if(similarity==3) lCn_1 = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, 1, 3, 1);
            if(similarity==4) lCn_1 = Similarity3_4(s1o, s2o, mu0, k0, v0, L0, 1, 4, 1);
            
            wei_empty(cc) = ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_empty.slice(cc).diag(0), Rinv_empty(0, 1, cc)) +
              lCn_1 + // Spatial part of similarity function
              log(Mdp) - log(CC);
            
          }
          
          arma::vec wei = arma::join_cols(wei_filled, wei_empty);
          
          // Trick to avoid all zero in the weights
          double max_wei = arma::max(wei);
          wei = wei - max_wei;
          wei = arma::trunc_exp(wei)/arma::sum(arma::trunc_exp(wei));
          
          // Sample newci based on wei (probability vector)
          int newci = sample_cpp(k_current + CC, 1, wei, true); // controlla l'output //
          
          // If newci is less than or equal to k_current, update existing cluster
          if (newci <= k_current) {
            // Point 2a: Update clu_current and nj_current for existing clusters
            clu_current(i) = newci;
            nj_current(newci - 1) += 1;
          } else {
            // Point 2b: Handle the case for new clusters
            int id_empty = newci - k_current;
            
            // Update the cluster for the new cluster
            clu_current(i) = k_current + 1;
            k_current += 1;
            
            nj_current.resize(k_current); // Resize nj_current to accommodate the new cluster
            nj_current(k_current - 1) = 1;
            
            phi_star.resize(k_current);
            phi_star(k_current - 1) = phi_empty(id_empty - 1);
            
            tau2_star.resize(k_current);
            tau2_star(k_current - 1) = tau2_empty(id_empty - 1);
            
            Phiinv_star = arma::join_slices(Phiinv_star, Phiinv_empty.slice(id_empty - 1));
            Rinv_star = arma::join_slices(Rinv_star, Rinv_empty.slice(id_empty - 1));
            
            phi_empty(id_empty - 1) = Rcpp::rbeta(1, a_phi, b_phi)[0];
            tau2_empty(id_empty - 1) = 1.0 / arma::randg<double>(arma::distr_param(a_tau, 1.0 / b_tau));
            
            Phiinv_empty.slice(id_empty - 1).diag(0).fill(1 + phi_empty(id_empty - 1) * phi_empty(id_empty - 1));
            Phiinv_empty.slice(id_empty - 1).diag(1).fill(-phi_empty(id_empty - 1));
            Phiinv_empty.slice(id_empty - 1).diag(-1).fill(-phi_empty(id_empty - 1));
            Phiinv_empty(0, 0, id_empty - 1) = 1.0;
            Phiinv_empty(T - 1, T - 1, id_empty - 1) = 1.0;
            
            Rinv_empty.slice(id_empty - 1) = 1.0 / tau2_empty(id_empty - 1) * Phiinv_empty.slice(id_empty - 1);
            
          }
          
        }
        
        phi(i) = phi_star(clu_current(i) - 1);
        tau2(i) = tau2_star(clu_current(i) - 1);
        
        Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv_star.slice(clu_current(i) - 1), arma::eye(T, T));
        Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv_star.slice(clu_current(i) - 1);
        
      }
      
    }
    
    ////////////////////
    // update phistar //
    ////////////////////
    
    arma::vec trans_phi_star = arma::trunc_log(phi_star) - arma::log1p(-phi_star);
    
    arma::vec trans_phi_star_prop = arma::mvnrnd(trans_phi_star, eta_phi * arma::eye(k_current, k_current), 1);
    
    arma::vec phi_star_prop = arma::trunc_exp(trans_phi_star_prop) / (1 + arma::trunc_exp(trans_phi_star_prop));
    
    arma::cube Phiinv_star_prop(T, T, k_current, arma::fill::zeros);
    arma::cube Rinv_star_prop(T, T, k_current, arma::fill::zeros);
    
    arma::vec log_target = (a_phi - 1) * arma::trunc_log(phi_star) +
      (b_phi - 1) * arma::log1p(-phi_star) +
      (trans_phi_star - 2 * arma::log1p(arma::trunc_exp(trans_phi_star)));
    
    arma::vec log_target_prop = (a_phi - 1) * arma::trunc_log(phi_star_prop) +
      (b_phi - 1) * arma::log1p(-phi_star_prop) +
      (trans_phi_star_prop - 2 * arma::log1p(arma::trunc_exp(trans_phi_star_prop)));
    
    for(int k = 0; k < k_current; k++){
      
      Phiinv_star_prop.slice(k).diag(0).fill(1 + phi_star_prop(k) * phi_star_prop(k));
      Phiinv_star_prop.slice(k).diag(1).fill(-phi_star_prop(k));
      Phiinv_star_prop.slice(k).diag(-1).fill(-phi_star_prop(k));
      Phiinv_star_prop(0, 0, k) = 1.0;
      Phiinv_star_prop(T - 1, T - 1, k) = 1.0;
      
      Rinv_star_prop.slice(k) = 1.0 / tau2_star(k) * Phiinv_star_prop.slice(k);
      
      for(int i = 0; i < n; i++){
        if(clu_current(i) - 1 == k){
          log_target_prop(k) += ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_star_prop.slice(k).diag(0), Rinv_star_prop(0, 1, k));
          log_target(k) += ldmvnorm_tridiag_equal_offdiag(w.col(i), arma::zeros(T), Rinv_star.slice(k).diag(0), Rinv_star(0, 1, k));
        }
      }
      
      double log_ratio = log_target_prop(k) - log_target(k);
      
      double log_alpha = std::min(0.0, log_ratio);
      
      double u = arma::randu();
      
      if(log(u) < log_alpha){
        
        acc_phi += 1;
        phi_star(k) = phi_star_prop(k);
        Phiinv_star.slice(k) = Phiinv_star_prop.slice(k);
        Rinv_star.slice(k) = Rinv_star_prop.slice(k);
        
      }
      
    }
    
    for(int i = 0; i < n; i++){
      phi(i) = phi_star(clu_current(i) - 1);
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv_star.slice(clu_current(i) - 1), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv_star.slice(clu_current(i) - 1);
    }
    
    /////////////////////
    // update tau2star //
    /////////////////////
    
    for(int k = 0; k < k_current; k++){
      double shape = a_tau + (nj_current(k) * T) / 2.0;
      double scale = b_tau;
      
      for(int i = 0; i < n; i++){
        if(clu_current(i) - 1 == k){
          scale += quad_tridiag(w.col(i), Phiinv_star.slice(k).diag(0), Phiinv_star(0, 1, k)) / 2.0;
        }
      }
      
      tau2_star(k) = 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / scale));
      
      Rinv_star.slice(k) = 1.0 / tau2_star(k) * Phiinv_star.slice(k);
      
    }
    
    for(int i = 0; i < n; i++){
      tau2(i) = tau2_star(clu_current(i) - 1);
      
      Q_n.slice(i) = invSigma_n.slice(i) - arma::solve(Sigma_n.slice(i) + sigma2(i) * sigma2(i) * Rinv_star.slice(clu_current(i) - 1), arma::eye(T, T));
      Vwinv_n.slice(i) = invSigma_n.slice(i) + Rinv_star.slice(clu_current(i) - 1);
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
      clustmemory.row(iter - burn) = clu_current.t();
      phistarsample[iter - burn] = phi_star;
      tau2starsample[iter - burn] = tau2_star;
      nj[iter - burn] = nj_current;
      k_clust(iter - burn) = k_current;
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
    Rcpp::Named("cluster") = clustmemory, 
    Rcpp::Named("nj") = nj,
    Rcpp::Named("phistar") = phistarsample,
    Rcpp::Named("accphi") = acc_phi/(double)runs,
    Rcpp::Named("tau2star") = tau2starsample,
    Rcpp::Named("k_clust") = k_clust,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("rem_obs") = (cts + 1)
  );
  
}