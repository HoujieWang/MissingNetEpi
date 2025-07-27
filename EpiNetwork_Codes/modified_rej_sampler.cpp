#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

#include <random>
#include <vector>
// We use the Rcpp namespace for convenience
using namespace Rcpp;

/**
 * survival_piecewise_const_exp_cpp
 *
 * Computes the survival function S(t) under a piecewise-constant hazard,
 * where 'intvl' contains sorted time boundaries (including 0 and possibly Inf),
 * and 'rate' contains the hazard rates for each interval.
 *
 * For t in the interval [ intvl[idx], intvl[idx+1] ),
 * the integrated hazard up to t is breaks[idx] + (t - intvl[idx]) * rate[idx].
 *
 * @param t time at which the survival is evaluated
 * @param intvl numeric vector of breakpoints (sorted ascending), e.g. c(0, b1, b2, ..., Inf)
 * @param rate  numeric vector of hazard rates (same length as number of intervals)
 * @return S(t) = exp( -integratedHazard(t) )
 *
 * Example: if intvl = (0, 2, 5, +Inf),
 *          and rate  = (0.1, 0.2, 0.3),
 *   - from 0 to 2  => hazard = 0.1
 *   - from 2 to 5  => hazard = 0.2
 *   - from 5+      => hazard = 0.3
 */
// [[Rcpp::export]]
double survival_piecewise_const_exp_cpp(double t,
                                        NumericVector intvl,
                                        NumericVector rate) {
  // Find which interval 't' falls in:
  // idx satisfies intvl[idx] <= t < intvl[idx+1]
  // We use std::lower_bound, then subtract 1 from the index
  // (be careful if t < intvl[0], you might get idx = -1)
  int idx = std::lower_bound(intvl.begin(), intvl.end(), t) - intvl.begin() - 1;
  
  // If t is smaller than intvl[0], or if idx < 0, survival is 1 by definition
  if (idx < 0) {
    return 1.0;
  }
  
  // Precompute the integrated hazard at each breakpoint in 'breaks'
  // breaks[i] = \int_{0}^{intvl[i]} hazard(s) ds
  NumericVector breaks(rate.size() + 1, 0.0);
  for (int i = 1; i < breaks.size(); ++i) {
    double dt = intvl[i] - intvl[i - 1];
    breaks[i] = breaks[i - 1] + rate[i - 1] * dt;
  }
  
  // If t >= intvl[last-1], idx could point to the last interval => clamp it
  if (idx >= (int)rate.size()) {
    idx = rate.size() - 1;
  }
  
  // Then total integrated hazard up to t is:
  //   breaks[idx] + (t - intvl[idx]) * rate[idx]
  double integratedHaz = breaks[idx] + (t - intvl[idx]) * rate[idx];
  return std::exp(-integratedHaz);
}

/**
 * trunc_piecewise_quantile_cpp
 *
 * R version was:
 *
 *   trunc_piecewise_quantile = function(u, breakpoints, lambda, T_max) {
 *     intervals = c(0, breakpoints, Inf)
 *     lhs = -log(1 - u*(1 - survival_piecewise_const_exp_cpp(T_max, intervals, lambda)))
 *     intvl = c(0, cumsum(diff(c(0, breakpoints)) * lambda[-length(lambda)]), Inf)
 *     intvl_idx = as.integer(cut(lhs, intvl, labels = 1:(length(breakpoints)+1)))
 *     (lhs - intvl[intvl_idx]) / lambda[intvl_idx] + c(0, breakpoints)[intvl_idx]
 *   }
 *
 * The logic here:
 * 1. intervals = (0, breakpoints..., Inf)
 * 2. lhs = -log( 1 - u * (1 - S(T_max)) )
 * 3. intvl = (0, cumsum(...), Inf)
 * 4. Find index such that intvl[idx] < lhs <= intvl[idx+1]
 * 5. Output = (lhs - intvl[idx]) / lambda[idx] + c(0, breakpoints)[idx]
 *
 * We handle the possibility that breakpoints might be empty. Then we have
 * just intervals = (0, +Inf), and a single hazard rate in lambda.
 */
// [[Rcpp::export]]
double trunc_piecewise_quantile_cpp(double u,
                                    NumericVector breakpoints,
                                    NumericVector lambda,
                                    double T_max) {
  // 1) Build 'intervals' for survival calculation => c(0, breakpoints, Inf)
  int K = breakpoints.size();
  NumericVector intervals(K + 2);
  intervals[0] = 0.0;
  for (int i = 0; i < K; i++) {
    intervals[i + 1] = breakpoints[i];
  }
  intervals[K + 1] = std::numeric_limits<double>::infinity();
  
  // 2) Compute S(T_max), then LHS = -log(1 - u*(1 - S(T_max)))
  double S_Tmax = survival_piecewise_const_exp_cpp(T_max, intervals, lambda);
  double lhs = -std::log(1.0 - u * (1.0 - S_Tmax));
  
  // If lhs <= 0, the quantile is 0 or near 0
  if (lhs <= 0.0) {
    return 0.0;
  }
  
  // 3) intvl = (0, cumsum(diff(c(0, breakpoints)) * lambda[-last]), Inf)
  //    Also shift = c(0, breakpoints).
  //    We'll build them step by step.
  
  // shift = c(0, breakpoints)
  NumericVector shift(K + 1);
  shift[0] = 0.0;
  for (int i = 0; i < K; i++) {
    shift[i + 1] = breakpoints[i];
  }
  
  // intvl
  // We'll accumulate in csum the integrated hazard for intervals up to i.
  NumericVector intvl(K + 2);
  intvl[0] = 0.0;
  
  double csum = 0.0;
  // We'll use lambda[i] for i up to lambda.size()-1 in the cumsum
  // because the last hazard extends to +Inf
  for (int i = 0; i < (int)lambda.size() - 1; i++) {
    double dt = shift[i + 1] - shift[i]; // breakpoints[i] - breakpoints[i-1]
    csum += dt * lambda[i];
    intvl[i + 1] = csum;
  }
  // The last cell is Inf
  intvl[K + 1] = std::numeric_limits<double>::infinity();
  
  // 4) Find idx such that intvl[idx] < lhs <= intvl[idx+1].
  // We'll just loop, but you could also use std::upper_bound / lower_bound with care.
  int idx = 0;
  for (int i = 0; i < (K + 1); i++) {
    if (lhs > intvl[i] && lhs <= intvl[i + 1]) {
      idx = i;
      break;
    }
  }
  
  // 5) Return (lhs - intvl[idx]) / lambda[idx] + shift[idx]
  double out = (lhs - intvl[idx]) / lambda[idx] + shift[idx];
  return out;
}

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

// You already have these two somewhere:
//  1) double trunc_piecewise_quantile_cpp(...);
//  2) DataFrame simulate_ctmc_cpp(...);

/**
 * simulate_ctmc3_cpp
 *
 * Mirrors the logic of your R function simulate_ctmc3().
 * It either calls simulate_ctmc_cpp() directly if first_state == last_state,
 * or it generates the initial transition time st_time, merges that single
 * event with subsequent events from simulate_ctmc_cpp(), and returns a
 * combined DataFrame of transitions.
 *
 * @param conn_rate    NumericVector of "connection" hazard rates (e.g., state=0 -> 1).
 * @param disconn_rate NumericVector of "disconnection" hazard rates (state=1 -> 0).
 * @param breakpoints  NumericVector of breakpoints (no need to include 0 or Inf).
 * @param T_max        Total simulation time (a double).
 * @param first_state  Starting state (0 or 1).
 * @param last_state   Desired final state (0 or 1).
 *
 * @return An R DataFrame with columns "time" and "event". The first row is the
 *         initial transition (if it occurs), and subsequent rows come from
 *         simulate_ctmc_cpp().
 */
// [[Rcpp::export]]
DataFrame simulate_ctmc3_cpp(NumericVector conn_rate,
                             NumericVector disconn_rate,
                             NumericVector breakpoints,
                             double T_max,
                             int first_state,
                             int last_state) 
{
  // -------------------------------------------------------------------------
  // 1) If first_state == last_state, just call simulate_ctmc_cpp() and return
  // -------------------------------------------------------------------------
  if (first_state == last_state) {
    // Return the result of simulate_ctmc_cpp() directly
    DataFrame out = simulate_ctmc_cpp(conn_rate, disconn_rate, breakpoints,
                                      T_max, first_state, last_state);
    return out;
  }
  
  // -------------------------------------------------------------------------
  // 2) Otherwise, we need to sample the first transition time st_time
  // -------------------------------------------------------------------------
  // 2.1) Pick the relevant piecewise hazard "lambda" based on first_state
  NumericVector lambda = (first_state == 0) ? conn_rate : disconn_rate;
  
  // 2.2) Generate a uniform random number 'u'
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  double u = runif(gen);
  
  // 2.3) If breakpoints is empty => pass an empty NumericVector to
  //      trunc_piecewise_quantile_cpp(). Otherwise, pass breakpoints as is.
  double st_time = 0.0;
  if (breakpoints.size() == 0) {
    NumericVector empty_bp(0); // an empty vector
    st_time = trunc_piecewise_quantile_cpp(u, empty_bp, lambda, T_max);
  } else {
    st_time = trunc_piecewise_quantile_cpp(u, breakpoints, lambda, T_max);
  }
  
  // -------------------------------------------------------------------------
  // 3) We now have the first transition time st_time. The state flips once.
  //    So the new current_time = st_time, new current_state = 1 - first_state.
  // -------------------------------------------------------------------------
  double current_time  = st_time;
  int    current_state = 1 - first_state;
  
  // -------------------------------------------------------------------------
  // 4) We shift and prune breakpoints to only those > st_time.
  //    Then subtract st_time so that the next sub-interval starts at 0.
  // -------------------------------------------------------------------------
  std::vector<double> new_breakpoints;
  new_breakpoints.reserve(breakpoints.size());
  for (double bp : breakpoints) {
    if (bp > st_time) {
      new_breakpoints.push_back(bp - st_time);
    }
  }
  
  // Turn new_breakpoints into an Rcpp NumericVector
  NumericVector shifted_bp = wrap(new_breakpoints);
  
  // -------------------------------------------------------------------------
  // 5) Subselect the portion of conn_rate / disconn_rate that matches
  //    the newly shortened breakpoints, just like R code:
  //
  //    nintvl = length(conn_rate)
  //    events = simulate_ctmc_cpp(
  //      conn_rate[(nintvl - length(shifted_bp)): nintvl],
  //      disconn_rate[(nintvl - length(shifted_bp)): nintvl],
  //      shifted_bp, T_max-st_time, current_state, last_state
  //    )
  // -------------------------------------------------------------------------
  int nintvl = conn_rate.size(); 
  int new_len = shifted_bp.size(); // how many breakpoints remain
  
  // We'll extract the last `new_len + 1` rates from each vector:
  // because if there are K breakpoints, we typically have K+1 rates.
  // So if new_len == 0 => we want the last 1 element, etc.
  int start_idx = nintvl - (new_len + 1);
  if (start_idx < 0) {
    start_idx = 0; // just in case
  }
  
  // Build sub-vectors
  NumericVector sub_conn_rate(new_len + 1), sub_disconn_rate(new_len + 1);
  for (int i = 0; i < (new_len + 1); i++) {
    sub_conn_rate[i]    = conn_rate[start_idx + i];
    sub_disconn_rate[i] = disconn_rate[start_idx + i];
  }
  
  // -------------------------------------------------------------------------
  // 6) Call simulate_ctmc_cpp(...) for the remainder of the time
  // -------------------------------------------------------------------------
  DataFrame events = simulate_ctmc_cpp(
    sub_conn_rate,
    sub_disconn_rate,
    shifted_bp,
    T_max - st_time,      // we only have (T_max - st_time) time left
    current_state,
    last_state
  );
  
  // -------------------------------------------------------------------------
  // 7) If we got > 0 rows back, shift the times by st_time
  // -------------------------------------------------------------------------
  int nrows = events.nrows(); 
  if (nrows > 0) {
    // Add st_time to each time in 'events'
    NumericVector e_time = events["time"];
    for (int i = 0; i < nrows; i++) {
      e_time[i] += st_time;
    }
    // Put the times back
    events["time"] = e_time;
  }
  
  // -------------------------------------------------------------------------
  // 8) We want to "rbind" the single initial transition with the subsequent
  //    events. The initial transition is (st_time, event = 1 - first_state).
  // -------------------------------------------------------------------------
  // We'll build two std::vectors for time and event, then combine.
  std::vector<double> final_times;
  std::vector<int>    final_events;
  
  // The first row: time = st_time, event = 1 - first_state
  final_times.push_back(st_time);
  final_events.push_back(1 - first_state);
  
  // Then append the rows from 'events'
  if (nrows > 0) {
    NumericVector e_time  = events["time"];
    IntegerVector e_event = events["event"];
    for (int i = 0; i < nrows; i++) {
      final_times.push_back(e_time[i]);
      final_events.push_back(e_event[i]);
    }
  }
  
  // -------------------------------------------------------------------------
  // 9) Construct and return the final DataFrame
  // -------------------------------------------------------------------------
  return DataFrame::create(
    Named("time")  = final_times,
    Named("event") = final_events
  );
}
// [[Rcpp::export]]
DataFrame run_simulations_cpp(List list_of_lists) {
  
  // We'll accumulate results in these vectors:
  std::vector<double> all_time;
  std::vector<int>    all_event;
  std::vector<int>    all_per1;
  std::vector<int>    all_per2;
  
  // Loop over each sub-list in the outer list
  for (int i = 0; i < list_of_lists.size(); i++) {
    // Extract the i-th sub-list
    List sublist = list_of_lists[i];
    
    // Each sub-list has 8 elements in a fixed format:
    //   1) conn_rate (NumericVector)
    //   2) disconn_rate (NumericVector)
    //   3) breakpoints (NumericVector)
    //   4) T_max (double)
    //   5) first_state (int)
    //   6) last_state (int)
    //   7) time_shift (double)
    //   8) per_vec (IntegerVector of length 2, for per1 and per2)
    
    NumericVector conn_rate    = sublist[0];
    NumericVector disconn_rate = sublist[1];
    NumericVector breakpoints  = sublist[2];
    double        T_max        = sublist[3];
    int           first_state  = sublist[4];
    int           last_state   = sublist[5];
    double        time_shift   = sublist[6];
    IntegerVector per_vec      = sublist[7]; // length 2 => [per1, per2]
    
    // Safety check (optional): make sure per_vec has length 2
    if (per_vec.size() != 2) {
      stop("per_vec must have length 2 in each sub-list!");
    }
    
    // 1) Call your existing simulate_ctmc3_cpp
    DataFrame sim_result = simulate_ctmc3_cpp(
      conn_rate,
      disconn_rate,
      breakpoints,
      T_max,
      first_state,
      last_state
    );
    
    // 2) Extract columns from sim_result
    NumericVector times  = sim_result["time"];
    IntegerVector events = sim_result["event"];
    
    // 3) For each row in sim_result, add to the accumulator vectors
    //    * times[r] + time_shift
    //    * events[r]
    //    * per_vec[0], per_vec[1]
    for (int r = 0; r < times.size(); r++) {
      all_time.push_back(times[r] + time_shift);
      all_event.push_back(events[r]);
      all_per1.push_back(per_vec[0]);
      all_per2.push_back(per_vec[1]);
    }
  }
  
  // Build one final data frame
  DataFrame out = DataFrame::create(
    Named("time")  = wrap(all_time),
    Named("event") = wrap(all_event),
    Named("per1")  = wrap(all_per1),
    Named("per2")  = wrap(all_per2)
  );
  
  return out;
}
// [[Rcpp::export]]
DataFrame run_simulations_outer_cpp(List link_info) {
  
  // Vectors to accumulate results
  std::vector<double> all_time;
  std::vector<int>    all_event;
  std::vector<int>    all_per1;
  std::vector<int>    all_per2;
  
  // Loop over each element of link_info (equivalent to lapply in R)
  for (int i = 0; i < link_info.size(); i++) {
    // Extract the i-th sublist
    List x = link_info[i];
    
    // Initialize a DataFrame to store the result
    DataFrame events;
    
    // Retry until the DataFrame has at least one row
    while (true) {
      events = run_simulations_cpp(x); // Call the existing C++ function
      
      // Check if the result has rows
      if (events.nrows() > 0) {
        break;
      }
    }
    
    // Extract columns from the events DataFrame
    NumericVector times  = events["time"];
    IntegerVector evts   = events["event"];
    IntegerVector per1   = events["per1"];
    IntegerVector per2   = events["per2"];
    
    // Append the results to the accumulator vectors
    for (int r = 0; r < times.size(); r++) {
      all_time.push_back(times[r]);
      all_event.push_back(evts[r]);
      all_per1.push_back(per1[r]);
      all_per2.push_back(per2[r]);
    }
  }
  
  // Combine all accumulated results into a single DataFrame
  DataFrame out = DataFrame::create(
    Named("time")  = wrap(all_time),
    Named("event") = wrap(all_event),
    Named("per1")  = wrap(all_per1),
    Named("per2")  = wrap(all_per2)
  );
  
  return out;
}