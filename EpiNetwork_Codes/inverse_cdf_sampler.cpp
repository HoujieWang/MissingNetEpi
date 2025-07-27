#include <Rcpp.h>
#include <cmath>
#include <functional>
using namespace Rcpp;

double bisection(std::function<double(double)> f, double lower, double upper, double tol = 1e-8, int max_iter = 100) {
  double a = lower, b = upper;
  double fa = f(a), fb = f(b);
  
  // Ensure the root is bracketed
  if (fa * fb > 0) {
    stop("The function must have different signs at the interval endpoints.");
  }
  
  double c, fc;
  for (int i = 0; i < max_iter; ++i) {
    c = (a + b) / 2.0; // Midpoint
    fc = f(c);
    
    if (std::abs(fc) < tol || (b - a) / 2.0 < tol) {
      return c; // Root found
    }
    
    // Update interval
    if (fa * fc < 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }
  
  stop("Bisection method did not converge within the maximum number of iterations.");
}

// Define the P function
NumericMatrix P(double t, double alpha, double omega) {
  NumericMatrix mat(2, 2);
  double exp_term = std::exp(-(alpha + omega) * t);
  
  mat(0, 0) = alpha + omega * exp_term;
  mat(0, 1) = omega - omega * exp_term;
  mat(1, 0) = alpha - alpha * exp_term;
  mat(1, 1) = omega + alpha * exp_term;
  
  return mat / (omega + alpha);
}

// Define the P_full function
NumericMatrix P_full(NumericVector times, NumericVector alpha, NumericVector omega) {
  NumericVector time_diff = diff(times);
  int n = time_diff.size();
  NumericMatrix P1 = P(time_diff[0], alpha[0], omega[0]);
  
  for (int k = 1; k < n; ++k) {
    NumericMatrix P2 = P(time_diff[k], alpha[k], omega[k]);
    NumericMatrix P_temp(2, 2);
    
    P_temp(0, 0) = P1(0, 0) * P2(0, 0) + P1(0, 1) * P2(1, 0);
    P_temp(0, 1) = P1(0, 0) * P2(0, 1) + P1(0, 1) * P2(1, 1);
    P_temp(1, 0) = P1(1, 0) * P2(0, 0) + P1(1, 1) * P2(1, 0);
    P_temp(1, 1) = P1(1, 0) * P2(0, 1) + P1(1, 1) * P2(1, 1);
    
    P1 = P_temp;
  }
  return P1;
}

// Define the survival_piecewise_const_exp_cpp function
// [[Rcpp::export]]
double survival_piecewise_const_exp_cpp(double t, NumericVector intvl, NumericVector rate) {
  int idx = std::lower_bound(intvl.begin(), intvl.end(), t) - intvl.begin() - 1;

  NumericVector breaks(rate.size() + 1, 0.0);
  for (int i = 1; i < breaks.size(); ++i) {
    breaks[i] = breaks[i - 1] + rate[i - 1] * (intvl[i] - intvl[i - 1]);
  }

  return std::exp(-(breaks[idx] + (t - intvl[idx]) * rate[idx]));
}

// Define the cdf function
// [[Rcpp::export]]
double cdf_cpp_diff(double t, double q, NumericVector full_intvl, NumericVector alpha, NumericVector omega, int net0, int net1) {
  // Determine the interval index for t
  int idx = std::distance(full_intvl.begin(), std::upper_bound(full_intvl.begin(), full_intvl.end(), t)) - 1;
  
  // Handle edge case: if t == upper, set idx to the last valid index
  if (idx >= full_intvl.size() - 1) {
    idx = full_intvl.size() - 2;  // Avoid out-of-range slicing
  }
  
  // Construct sub_intvl
  NumericVector sub_intvl = NumericVector::create(t);
  if (idx + 1 < full_intvl.size()) {
    NumericVector tail = full_intvl[Range(idx + 1, full_intvl.size() - 1)];
    sub_intvl = NumericVector(sub_intvl.size() + tail.size());
    sub_intvl[0] = t;
    std::copy(tail.begin(), tail.end(), sub_intvl.begin() + 1);
  }
  
  // Debugging: Print sub_intvl
  // Rcpp::Rcout << "sub_intvl: ";
  // for (double v : sub_intvl) Rcpp::Rcout << v << " ";
  // Rcpp::Rcout << std::endl;
  
  // Compute P_full for sub_intvl and full_intvl
  NumericMatrix P_sub = P_full(sub_intvl, alpha[Range(idx, alpha.size() - 1)], omega[Range(idx, omega.size() - 1)]);
  NumericMatrix P_full_matrix = P_full(full_intvl, alpha, omega);
  
  // Compute survival_piecewise_const_exp_cpp
  NumericVector truncated_full_intvl = full_intvl[Range(0, full_intvl.size() - 2)];
  truncated_full_intvl.push_back(R_PosInf);
  
  NumericVector param;
  if (net0 == 0){ // 0 means connection so next event is disconnection
    param = omega;
  }else{
    param = alpha;
  }
  double survival = survival_piecewise_const_exp_cpp(t, truncated_full_intvl, param);
  
  // Calculate result
  double result = 1 - (P_sub(net0, net1) / P_full_matrix(net0, net1)) * survival - q;
  
  return result;
}
// [[Rcpp::export]]
double cdf_cpp_same(double t, double q, NumericVector full_intvl, NumericVector alpha, NumericVector omega, int net0, int net1) {
  // Determine the interval index for t
  int idx = std::distance(full_intvl.begin(), std::upper_bound(full_intvl.begin(), full_intvl.end(), t)) - 1;
  
  // Handle edge case: if t == upper, set idx to the last valid index
  if (idx >= full_intvl.size() - 1) {
    idx = full_intvl.size() - 2;  // Avoid out-of-range slicing
  }
  
  // Construct sub_intvl: c(t, full_intvl[(idx+1):length(full_intvl)])
  NumericVector sub_intvl = NumericVector::create(t);
  if (idx + 1 < full_intvl.size()) {
    NumericVector tail = full_intvl[Range(idx + 1, full_intvl.size() - 1)];
    sub_intvl = NumericVector(sub_intvl.size() + tail.size());
    sub_intvl[0] = t;
    std::copy(tail.begin(), tail.end(), sub_intvl.begin() + 1);
  }
  
  // Debug: Print sub_intvl
  // Rcpp::Rcout << "sub_intvl: ";
  // for (double v : sub_intvl) Rcpp::Rcout << v << " ";
  // Rcpp::Rcout << std::endl;
  
  // Select the parameter based on net0
  // NumericVector param;
  // if (net0 == 0){
  //   param = alpha;
  // }else{
  //   param = omega;
  // }
  // 
  // Compute Tmax
  double Tmax = full_intvl[full_intvl.size() - 1];
  
  // Compute p_above
  NumericVector truncated_full_intvl = full_intvl[Range(0, full_intvl.size() - 2)];
  truncated_full_intvl.push_back(R_PosInf);
  
  
  double survival_Tmax;
  double survival_t;
  if (net0 == 0){
     survival_Tmax = survival_piecewise_const_exp_cpp(Tmax, truncated_full_intvl, omega);
     survival_t = survival_piecewise_const_exp_cpp(t, truncated_full_intvl, omega);
  }else{
     survival_Tmax = survival_piecewise_const_exp_cpp(Tmax, truncated_full_intvl, alpha);
     survival_t = survival_piecewise_const_exp_cpp(t, truncated_full_intvl, alpha);
  }
  double p_above = survival_Tmax / P_full(full_intvl, alpha, omega)(net0, net0);
  // Debug: Print p_above components
  // Rcpp::Rcout << "Tmax: " << Tmax << ", survival_Tmax: " << survival_Tmax << ", truncated_full_intvl: " << truncated_full_intvl
  //             << ", P_full_matrix(net0, net0): " << P_full(full_intvl, alpha, omega)(net0, net0)
  //             << ", p_above: " << p_above << std::endl;
  
  // Compute P_full for sub_intvl
  NumericMatrix P_sub = P_full(sub_intvl, alpha[Range(idx, alpha.size() - 1)], omega[Range(idx, omega.size() - 1)]);
  NumericMatrix P_full_matrix = P_full(full_intvl, alpha, omega);
  
  // Debug: Print P_sub and P_full_matrix
  // Rcpp::Rcout << "P_sub(net0, net0): " << P_sub(net0, net0) << ", P_full_matrix(net0, net0): " << P_full_matrix(net0, net0) << std::endl;
  // 
  
  
  // Debug: Print survival_t
  // Rcpp::Rcout << "survival_t: " << survival_t << std::endl;
  
  // Calculate numerator
  double numerator = (P_sub(net0, net0) / P_full_matrix(net0, net0)) * survival_t - p_above;
  
  // Debug: Print numerator
  // Rcpp::Rcout << "numerator: " << numerator << std::endl;
  
  // Calculate result
  double result = 1 - (numerator / (1 - p_above)) - q;
  
  // Debug: Print final result components
  // Rcpp::Rcout << "result: " << result << ", 1 - numerator / (1 - p_above): "
  //             << 1 - (numerator / (1 - p_above)) << ", q: " << q << std::endl;
  // 
  return result;
}

int if_impute(double q, NumericVector full_intvl, NumericVector alpha, NumericVector omega, int net0, int net1) {

  // Compute Tmax
  double Tmax = full_intvl[full_intvl.size() - 1];
  
  // Compute p_above
  NumericVector truncated_full_intvl = full_intvl[Range(0, full_intvl.size() - 2)];
  truncated_full_intvl.push_back(R_PosInf);
  
  
  double survival_Tmax;
  if (net0 == 0){
    survival_Tmax = survival_piecewise_const_exp_cpp(Tmax, truncated_full_intvl, omega);
  }else{
    survival_Tmax = survival_piecewise_const_exp_cpp(Tmax, truncated_full_intvl, alpha);
  }
  int result = (q > survival_Tmax) ? 1 : 0;
  return result;
}


// Define the main loop
NumericVector sample_event_times(NumericMatrix intvls, NumericVector u, NumericMatrix alpha, NumericMatrix omega, 
                                 IntegerVector net0, IntegerVector net1, bool diff) {
  int n = intvls.nrow();
  NumericVector event_times_sampled(n);
  
  for (int i = 0; i < n; ++i) {
    // Rcpp::Rcout << "Iteration i: " << i << std::endl;
    // Extract the i-th row of alpha and omega as vectors
    NumericVector alpha_row = alpha(i, _);
    NumericVector omega_row = omega(i, _);
    
    double lower = intvls(i, 0);
    double upper = intvls(i, intvls.ncol() - 1);
    
    // Rcpp::Rcout << "Row " << i << ": lower = " << lower << ", upper = " << upper << std::endl;
    
    // Skip rows where lower >= upper
    // if (lower >= upper) {
    //   // Rcpp::Rcout << "Skipping invalid row " << i << ": lower >= upper" << std::endl;
    //   continue;
    // }
    
    // Define a lambda function for cdf that uses the current row of alpha and omega

    auto cdf_lambda1 = [&](double t) -> double {
      return cdf_cpp_diff(t, u[i], intvls(i, _), alpha_row, omega_row, net0[i], net1[i]);
    };
   
    auto cdf_lambda2 = [&](double t) -> double {
      return cdf_cpp_same(t, u[i], intvls(i, _), alpha_row, omega_row, net0[i], net1[i]);
    };
    
    double root;
    // Use a custom bisection method for root finding
    if (diff){
      root = bisection(cdf_lambda1, lower, upper);
    }else{
      root = bisection(cdf_lambda2, lower, upper);
    }
    event_times_sampled[i] = root;
  }
  
  return event_times_sampled;
}
// [[Rcpp::export]]
IntegerVector check_all_if_impute(NumericMatrix intvls, NumericVector u, NumericMatrix alpha, NumericMatrix omega, 
                                 IntegerVector net0, IntegerVector net1) {
  int n = intvls.nrow();
  IntegerVector event_times_sampled(n);
  
  for (int i = 0; i < n; ++i) {
    // Rcpp::Rcout << "Iteration i: " << i << std::endl;
    // Extract the i-th row of alpha and omega as vectors
    NumericVector alpha_row = alpha(i, _);
    NumericVector omega_row = omega(i, _);
    
    int out;
    out = if_impute(u[i], intvls(i, _), alpha_row, omega_row, net0[i], net1[i]);
  
    event_times_sampled[i] = out;
  }
  
  return event_times_sampled;
}

// [[Rcpp::export]]
NumericVector cpp_sample_event_times(NumericMatrix intvls, NumericVector u, NumericMatrix alpha, NumericMatrix omega, IntegerVector net0, IntegerVector net1, bool diff) {
  return sample_event_times(intvls, u, alpha, omega, net0, net1, diff);
}