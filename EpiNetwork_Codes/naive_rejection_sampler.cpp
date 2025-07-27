#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <random>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame simulate_ctmc_cpp(NumericVector conn_rate,
                            NumericVector disconn_rate,
                            NumericVector breakpoints,
                            double T_max,
                            int first_state,
                            int last_state) {
  
  // 1. Build the intervals vector: [0, breakpoints..., Inf]
  NumericVector intervals(1, 0.0);  // Start with 0
  for (double bp : breakpoints) {
    intervals.push_back(bp);
  }
  intervals.push_back(R_PosInf);    // Append +Inf at the end
  
  // 2. Prepare output containers
  std::vector<double> times;  // transition times
  std::vector<int>    states; // states at each transition
  
  // 3. Initialize simulation
  double current_time  = 0.0;
  int    current_state = first_state;
  
  // 4. Set up C++ random number generation (uniform [0,1])
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  
  // 5. Outer loop: repeat until final state is reached
  while (true) {
    
    // Inner loop: keep simulating until we exceed T_max
    while (current_time < T_max) {
      
      // Pick the current rate vector depending on the state
      NumericVector rate = (current_state == 0) ? conn_rate : disconn_rate;
      
      // --- A. Find which sub-interval `current_time` belongs to ---
      int cur_index = 0; 
      for (int i = 1; i < intervals.size(); ++i) {
        if (current_time >= intervals[i - 1] && current_time < intervals[i]) {
          cur_index = i - 1;
          break;
        }
      }
      
      // --- B. Construct intvl_current ---
      //     R code does: intvl_current = c(current_time, intervals[(cur_index+1):end]) - current_time
      //     The first entry is 0, then intervals[i] - current_time for i in [cur_index+1..end]
      int N = intervals.size();
      int n_intvl = N - cur_index;  // length of intvl_current
      NumericVector intvl_current(n_intvl);
      
      // First element: current_time - current_time = 0
      intvl_current[0] = 0.0;
      // Subsequent elements: intervals[i] - current_time
      for (int i = cur_index + 1; i < N; i++) {
        int index_local = i - (cur_index + 1) + 1; // offset in intvl_current
        intvl_current[index_local] = intervals[i] - current_time;
      }
      
      // --- C. Construct rate_current ---
      //     R code does: rate_current = rate[cur_index : length(rate)]
      int M = rate.size();
      int n_rate = M - cur_index;
      NumericVector rate_current(n_rate);
      for (int i = cur_index; i < M; i++) {
        rate_current[i - cur_index] = rate[i];
      }
      
      // --- D. Generate random exponential scale: q = -log(1 - runif(1)) ---
      double u = runif(gen);
      double q = -std::log(1.0 - u); // same as rexp(1) in R
      
      // --- E. Compute the cumulative 'breaks' from piecewise intervals
      //     breaks = cumsum( c(0, diff(intvl_current) * rate_current ) )
      //     Length = n_intvl because we have 1 + (n_intvl-1)
      if (n_intvl - 1 != n_rate) {
        // In some cases, your piecewise vectors might have different lengths.
        // The below logic assumes n_intvl-1 == n_rate. If not, handle it gracefully or adjust your input rates.
      }
      
      int n_breaks = n_intvl; 
      NumericVector breaks(n_breaks);
      breaks[0] = 0.0;
      for (int i = 1; i < n_breaks; i++) {
        double delta = intvl_current[i] - intvl_current[i - 1];   // diff(intvl_current)[i-1]
        double val   = delta * rate_current[i - 1];               // * rate_current
        breaks[i]    = breaks[i - 1] + val;
      }
      
      // --- F. Determine which sub-interval the random q falls into
      //     next_intvl_idx = as.integer(cut(q, breaks, labels=1:(length(breaks)-1)))
      int next_intvl_idx = n_breaks - 2; // default to last if not found
      for (int i = 1; i < n_breaks; i++) {
        if (q >= breaks[i - 1] && q < breaks[i]) {
          next_intvl_idx = i - 1;
          break;
        }
      }
      
      // --- G. Compute new_time
      //     new_time = 1/rate_current[next_intvl_idx] * (q - breaks[next_intvl_idx])
      //                + intvl_current[next_intvl_idx] + current_time
      double rate_val = rate_current[next_intvl_idx];
      double new_time = (1.0 / rate_val) * (q - breaks[next_intvl_idx]) 
        + intvl_current[next_intvl_idx] 
      + current_time;
      
      // --- H. Check if we exceed T_max
      if (new_time > T_max) {
        break;  // exit inner while loop
      }
      
      // --- I. Record the transition ---
      current_time  = new_time;
      current_state = 1 - current_state;  // toggle state 0<->1
      times.push_back(current_time);
      states.push_back(current_state);
    } // end while(current_time < T_max)
    
    // --- J. Check final-state condition
    if (current_state == last_state) {
      // We have reached the required final state, so stop
      break;
    } else {
      // Reset simulation if we didn't end in the required final state
      current_time  = 0.0;
      current_state = first_state;
      times.clear();
      states.clear();
    }
  } // end while(true)
  
  // Return a data frame with two columns: "time" and "event"
  return DataFrame::create(
    Named("time")  = wrap(times),
    Named("event") = wrap(states)
  );
}