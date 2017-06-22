// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List funcKM2_eigen(   Eigen::Map<Eigen::VectorXd> t_idx,
                      Eigen::Map<Eigen::VectorXd> status,
                      Eigen::Map<Eigen::VectorXd> x,
                      Eigen::Map<Eigen::VectorXd> weight
                   ){
  
  using namespace Eigen;
  typedef Array<bool,Dynamic,1> ArrayXb;
  
  // Initiate value
  int n = t_idx.size();
  VectorXd Y(n); Y.fill(0);  //n.risk
  VectorXd N(n); N.fill(0);  //n.event
  VectorXd C(n); C.fill(0);  //n.censor
  VectorXd S(n); S.fill(0);  //surv
  VectorXd H(n); H.fill(0);  //forSE
  VectorXd D(n); D.fill(0);  //forSE
  VectorXd E(n); E.fill(0);  //SE

  Y[0] = weight.sum();
  N[0] = 0;
  C[0] = 0;
  S[0] = 1;
  H[0] = 0;
  D[0] = 0;
  E[0] = 0;
  
  ArrayXb term0(n);
  ArrayXb term1(n);
  ArrayXb term2(n);
  
  for(int i = 1; i < n; i++){
    Y[i] = Y[i-1] - N[i-1] - C[i-1];
    
    term0 = x.array() == t_idx[i]; // x == t_idx[i]
    term1 = status.array() == 1; // status == 1
    term2 = status.array() == 0; // status == 0
    
    // N
    if( (term0 && term1).cast <double>().sum() > 0){
      
      // sum(x==t_idx[i] & status==1)* weight[which(x==t_idx[i] & status==1)[1]] 
      int j = 0; 
      while( ! (term0[j] && term1[j]) ){
        j++;
      }
    
      N[i] = (term0 && term1).cast <double>().sum() * weight[j];
      
    }else{
      
      N[i] = 0;
    }
    
    // C
    if( (term0 && term2).cast <double>().sum() > 0){
      
      // sum(x==t_idx[i] & status==0)* weight[which(x==t_idx[i] & status==0)[1]] 
      int j = 0; 
      while( ! (term0[j] && term2[j])){
        j++;
      }
      
      C[i] = (term0 && term2).cast <double>().sum() * weight[j];
      
    }else{
      
      C[i] = 0;
    }
    
    if(Y[i]<0){Y[i] = 0;}
    
    if(Y[i]==0){
      S[i] = S[i-1];
    }else{
      S[i] = S[i-1]*(1-(N[i]/Y[i]));
    }
    
    if(Y[i]*(Y[i]-N[i])==0){
      H[i] = 0;
    }else{
      H[i] = N[i]/(Y[i]*(Y[i]-N[i]));
    }
    
    if(S[i]<0) S[i] = 0;
    
    D[i] = H.segment(1,i).sum();
    E[i] = sqrt(S[i]*S[i]*D[i]);
    
  }
  
  

  
  // return 
  return Rcpp::List::create(Rcpp::Named("t_idx") = t_idx,
                            Rcpp::Named("n_risk") = Y,
                            Rcpp::Named("n_event") = N,
                            Rcpp::Named("n_censor") = C,
                            Rcpp::Named("surv") = S,
                            Rcpp::Named("SE") = E
  );
                              
} 

