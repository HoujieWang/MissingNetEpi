#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec big_sumer_cpp2(arma::vec event_time, 
                        arma::vec event_type,
                        arma::vec event_p1,
                        arma::vec event_p2,
                        arma::mat adjmat, 
                        arma::vec I0) {
  int N = adjmat.n_rows;
  arma::vec epid(N); 
  for (int i = 0; i < I0.n_elem; ++i){
    epid(I0(i)-1) = 1;
  }
  // epid[I0-1] = 1;
  arma::uvec susceptible = find(epid == 0); 
  arma::uvec infected = find(epid == 1); 
  arma::uvec recovered = find(epid == -1); 
  int N_S = susceptible.n_elem; 
  int N_I = infected.n_elem;
  int N_R = sum(event_type == 2);
  arma::vec Mmax(3);
  Mmax[0] = N_S*(N_S-1)/2; Mmax[1] = N_S * N_I; Mmax[2] =  N_I*(N_I-1)/2;
  arma::mat tmp_mat = adjmat(susceptible, susceptible);
  arma::vec M = {accu(adjmat(susceptible, susceptible))/2, 
                 accu(adjmat(susceptible, infected)), 
                 accu(adjmat(infected, infected))/2};
  int N_SI = M[1];
  // int n_E = sum(event_type == 1);
  // int n_R = sum(event_type == 2);
  
  double t_cur = 0;
  int z; int I_p1; int R_p1; int S_p1; int p1; int p2;
  double t_next;
  double del_t;
  arma::vec big_sums(9);
  arma::vec C(3);
  arma::vec D(3);
  int n_event = event_type.n_elem;
  arma::rowvec p1_row;
  for(int r = 0; (r < n_event); ++r) {
    t_next = event_time[r];
    del_t = t_next - t_cur;
    // update "big_sums"
    big_sums[0] = big_sums[0] + N_SI*del_t;
    big_sums[1] = big_sums[1] + N_I*del_t;
    big_sums[2] = big_sums[2] + (Mmax[0] - M[0])*del_t;
    big_sums[3] = big_sums[3] + (Mmax[1] - M[1])*del_t;
    big_sums[4] = big_sums[4] + (Mmax[2] - M[2])*del_t;
    big_sums[5] = big_sums[5] + M[0]*del_t;
    big_sums[6] = big_sums[6] + M[1]*del_t;
    big_sums[7] = big_sums[7] + M[2]*del_t;
    big_sums[8] = big_sums[8] +  N_S*del_t;
    
    z = event_type[r];
    if(z == 1 || z == 0){
      p1 = event_p1[r];
      // figure out neighborhood
      p1_row = adjmat.row(p1-1);
      susceptible = find(epid == 0); 
      infected = find(epid == 1); 
      recovered = find(epid == -1);
      S_p1 = accu(p1_row(susceptible));
      I_p1 = accu(p1_row(infected));
      R_p1 = accu(p1_row(recovered));
      
      // bookkeeping
      epid[p1-1] = 1;
      N_S -= 1;
      N_I += 1;
      
      //Mmax = c((N_S+N_R)*(N_S+N_R-1)/2, (N_S+N_R) * N_I, N_I*(N_I-1)/2)
      Mmax[0] = (N-N_I)*(N-N_I-1)/2;
      Mmax[1] = (N-N_I) * N_I;
      Mmax[2] = N_I*(N_I-1)/2;
      
      M[0] -= (S_p1 + R_p1);
      M[1] += (S_p1 + R_p1 - I_p1);
      M[2] += I_p1 ;
      N_SI += (S_p1 - I_p1);
    }else if(z == 2){
      p1 = event_p1[r];
      p1_row = adjmat.row(p1-1);
      susceptible = find(epid == 0); 
      infected = find(epid == 1); 
      recovered = find(epid == -1);
      S_p1 = accu(p1_row(susceptible));
      I_p1 = accu(p1_row(infected));
      R_p1 = accu(p1_row(recovered));
      
      // bookkeeping
      N_I -= 1;
      
      epid[p1-1] = -1;
      N_R += 1;
      M[0] += (S_p1 + R_p1);
      M[1] += (-S_p1 - R_p1 + I_p1);
      M[2] -= I_p1;
      
      Mmax[0] = (N-N_I)*(N-N_I-1)/2;
      Mmax[1] = (N-N_I) * N_I;
      Mmax[2] = N_I*(N_I-1)/2;
      N_SI -= S_p1;
    }else{
      p1 = event_p1[r];
      p2 = event_p2[r];
      if ((3 <= z) & (z <= 5)){
        adjmat(p1-1, p2-1) = 1; 
        adjmat(p2-1, p1-1) = 1;
        
        if ((epid[p1-1]!=1) & (epid[p2-1]!=1)){
          // S-S/S-R/R-R type
          C[0] += 1;
          M[0] += 1;
        }else if ((epid[p1-1]==1) & (epid[p2-1]==1)){
          // I-I type
          C[2] += 1;
          M[2] += 1;
        }else if ((epid[p1-1]!=-1) & (epid[p2-1]!=-1)){
          // S-I type: neither is recovered
          C[1] += 1;
          M[1] += 1;
          N_SI += 1;
        }else{
          // R-I type
          C[1] += 1;
          M[1] += 1;
        }
        
        
      }else{
        adjmat(p1-1, p2-1) = 0; 
        adjmat(p2-1, p1-1) = 0;
        // bookkeeping
        if ((epid[p1-1]!=1) & (epid[p2-1]!=1)){
          // S-S/S-R/R-R type
          D[0] += 1;
          M[0] -= 1;
        }else if ((epid[p1-1]==1) &( epid[p2-1]==1)){
          // I-I type
          D[2] += 1;
          M[2] -= 1;
        }else if ((epid[p1-1]!=-1) & (epid[p2-1]!=-1)){
          // S-I type
          D[1] += 1;
          M[1] -= 1;
          N_SI -= 1;
        }else{  
          // R-I type
          D[1] += 1;
          M[1] -= 1;
        }
      }
    }
    t_cur = t_next;
  }
  return arma::join_vert(big_sums, arma::join_vert(C, D));
  // return big_sums;
}
