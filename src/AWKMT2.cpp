
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

struct KM{
  Eigen::ArrayXd t_idx;  //ordered unique observed time
  Eigen::ArrayXd Y;  //n.risk
  Eigen::ArrayXd N;  //n.event
  Eigen::ArrayXd C;  //n.censor
  Eigen::ArrayXd S;  //surv
  Eigen::ArrayXd E;  //SE
};

struct V1V2{
  Eigen::ArrayXd V1;         // V1 test statistics
  Eigen::ArrayXd V2;         // V2 test statistics
  Eigen::ArrayXd surv_diff;  // survival difference
};


NumericVector sort_rcpp(NumericVector x) {
  std::vector<double> tmp = Rcpp::as< std::vector<double> > (x);
  std::sort(tmp.begin(), tmp.end());
  return wrap(tmp);
}

IntegerVector rank_rcpp(NumericVector x) {
  return match(x, sort_rcpp(x));
}

KM funcKM2_eigen(   Eigen::ArrayXd t_idx,
                    Eigen::ArrayXd status,
                    Eigen::ArrayXd x,
                    Eigen::ArrayXd weight
){

  using namespace Eigen;
  typedef Array<bool,Dynamic,1> ArrayXb;

  // Initiate value
  int n = t_idx.size();
  ArrayXd Y(n); Y.fill(0);  //n.risk
  ArrayXd N(n); N.fill(0);  //n.event
  ArrayXd C(n); C.fill(0);  //n.censor
  ArrayXd S(n); S.fill(0);  //surv
  ArrayXd H(n); H.fill(0);  //forSE
  ArrayXd D(n); D.fill(0);  //forSE
  ArrayXd E(n); E.fill(0);  //SE

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

    term0 = x == t_idx[i]; // x == t_idx[i]
    term1 = status == 1; // status == 1
    term2 = status == 0; // status == 0

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
  // return Rcpp::List::create(Rcpp::Named("t_idx") = t_idx,
                               //                           Rcpp::Named("n_risk") = Y,
                               //                           Rcpp::Named("n_event") = N,
                               //                           Rcpp::Named("n_censor") = C,
                               //                           Rcpp::Named("surv") = S,
                               //                           Rcpp::Named("SE") = E
                               // );
  KM res;
  res.t_idx = t_idx;
  res.Y = Y;
  res.N = N;
  res.C = C;
  res.S = S;
  res.E = E;
  return res;

}

V1V2 funcAWKMT2_eigen(   Eigen::ArrayXd t_idx,
                         Eigen::ArrayXd status0,
                         Eigen::ArrayXd status1,
                         Eigen::ArrayXd x0,
                         Eigen::ArrayXd x1,
                         double tau1,
                         double tau2,
                         Eigen::ArrayXd crange,
                         CharacterVector test,
                         CharacterVector type,
                         Eigen::ArrayXd obs_survdiff
){
  using namespace Eigen;
  typedef Array<bool,Dynamic,1> ArrayXb;

  int n  = t_idx.size();
  int n0 = x0.size();
  int n1 = x1.size();
  int n_crange = crange.size();

  ArrayXd wt0(n0); wt0.fill(1);
  ArrayXd wt1(n1); wt1.fill(1);

  if( type[0] == ("perturbation") ){
    wt1 = as<Map<ArrayXd> >( rexp(n1) );
    wt0 = as<Map<ArrayXd> >( rexp(n0) );
  }

  //-- Get stats1 using funcKM2 --
    KM wk0; wk0 = funcKM2_eigen(t_idx, status0, x0, wt0);
    KM wk1; wk1 = funcKM2_eigen(t_idx, status1, x1, wt1);

    //-- Get stats2 --
      ArrayXd dt_diff(n);
    dt_diff[0] = 0;
    for(int i = 1; i < n; i++){
      dt_diff[i] = t_idx[i] - t_idx[i-1];
    }
    ArrayXd dt_jump(n); dt_jump = wk0.N + wk1.N; // "n_event1" + "n_event0"

    ArrayXb tau_idx(n); tau_idx = (t_idx > tau1) && (t_idx <=tau2);
    ArrayXd survdiff(n); survdiff = wk1.S - wk0.S; // surv1 - surv0
    ArrayXd survdiff_se(n); survdiff_se =  (wk1.E.pow(2) + wk0.E.pow(2)).sqrt();

    if(type[0] == ("perturbation")){
        survdiff = survdiff - obs_survdiff;
    }

    // -- Get Z1(one-sided) and Z2(two-sided) --
      ArrayXd z1(n); z1 = (survdiff_se == 0).select(0, survdiff/survdiff_se);
    ArrayXd z2(n); z2 = z1.abs();

    ArrayXd V1(n_crange);
    ArrayXd V2(n_crange);


    if( test[0] == "1_side"){
      for(int i = 0; i < crange.size(); i++){
        V1[i] = ( z1.max(crange[i]) * z1 * dt_diff * tau_idx.cast<double>()  ).sum();
        V2[i] = ( z1.max(crange[i]) * z1 * dt_jump * tau_idx.cast<double>()  ).sum();
      }
    }

    if( test[0] == "2_side"){
      for(int i = 0; i < crange.size(); i++){
        V1[i] = ( z2.array().max(crange[i]) * z2 * dt_diff * tau_idx.cast<double>()  ).sum();
        V2[i] = ( z2.array().max(crange[i]) * z2 * dt_jump * tau_idx.cast<double>()  ).sum();
      }
    }

    V1V2 res;
    res.V1 = V1;
    res.V2 = V2;
    res.surv_diff = survdiff;
    return res;
}

// [[Rcpp::export]]
List AWKMT2_eigen(       Eigen::Map<Eigen::ArrayXd> t_idx,
                         Eigen::Map<Eigen::ArrayXd> status0,
                         Eigen::Map<Eigen::ArrayXd> status1,
                         Eigen::Map<Eigen::ArrayXd> x0,
                         Eigen::Map<Eigen::ArrayXd> x1,
                         double tau1,
                         double tau2,
                         Eigen::Map<Eigen::ArrayXd> crange,
                         CharacterVector test,
                         CharacterVector type,
                         int nmethod
){
  using namespace Eigen;

  int n_crange = crange.size();
  int n0 = x0.size();
  int n1 = x1.size();
  int n01; n01 = n0 + n1;
  NumericVector ind(n01);
  for(int i =1; i < n01; i++){ ind[i] = i;}
  NumericVector perm_index(n01);

  V1V2 obs_ours; obs_ours = funcAWKMT2_eigen(t_idx, status0, status1, x0, x1, tau1, tau2, crange, test, "observe", t_idx);
  V1V2 tmp;

  MatrixXd V1(nmethod, n_crange);
  MatrixXd V2(nmethod, n_crange);

  ArrayXd x0_new(n0);
  ArrayXd x1_new(n1);
  ArrayXd x(n01); x << x0, x1;

  ArrayXd status0_new(n0);
  ArrayXd status1_new(n1);
  ArrayXd status(n01); status << status0, status1;

  for(int i = 0; i < nmethod; i++){
    if(type[0] == ("permutation")){

      perm_index = rank_rcpp(runif(n01)) - 1;

      int j0 = 0;
      int j1 = 0;
      for(int j = 0; j < n01; j++){
        if( perm_index[j] < n0){
          x0_new[j0] = x[j];
          status0_new[j0] = status[j];
          j0++;
        }else{
          x1_new[j1] = x[j];
          status1_new[j1] = status[j];
          j1++;
        }
      }

      tmp = funcAWKMT2_eigen(t_idx, status0_new, status1_new, x0_new, x1_new, tau1, tau2, crange, test, "permutation", obs_ours.surv_diff);


      }

    if(type[0] == ("perturbation")){
      tmp = funcAWKMT2_eigen(t_idx, status0, status1, x0, x1, tau1, tau2, crange, test, "perturbation", obs_ours.surv_diff );
    }

    V1.row(i) = tmp.V1;
    V2.row(i) = tmp.V2;

  }



  return Rcpp::List::create(  Rcpp::Named("V1") = obs_ours.V1,
                              Rcpp::Named("V2") = obs_ours.V2,
                              Rcpp::Named("V1_per") = V1,
                              Rcpp::Named("V2_per") = V2

  );

}
