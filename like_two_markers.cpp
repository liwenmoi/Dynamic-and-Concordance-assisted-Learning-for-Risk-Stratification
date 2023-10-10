#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector fracpoly(double t){
  // must put 1 as the first component.
  NumericVector frac(6);
  t = t + 1;
  frac[0]=1; frac[1]=log(t); frac[2]=sqrt(t); frac[3]=1/sqrt(t); frac[4]=t; frac[5]=1/(t);
  
  return frac;
}


// [[Rcpp::export]]
IntegerVector gen_seq_C(int n1, int n2){
  // gives a sequence of number from n1 to n2
  IntegerVector Aseq(n2-n1+1);
  Aseq[0] = n1;
  for (int i=1; i<n2-n1+1; i++){
    Aseq[i] = n1 + i;
  }
  
  return Aseq;
}


IntegerVector which_C(IntegerVector x, int y){
  // gives the index starting from 0
  IntegerVector v = gen_seq_C(0, x.size()-1);
  IntegerVector idx = v[x == y];
  return idx;
}

// [[Rcpp::export]]
IntegerVector test_C(IntegerVector dat_id, IntegerVector dat_delta, IntegerVector id){
  IntegerVector id_ev = dat_id[dat_delta == 1];
  
  // step 2: for each subject who had event, 
  
  int id_i = id_ev[0];
  return ( which_C(id, id_i));         
  
}


//[[Rcpp::export]]
double KerSmooth(double s, double t, double h_n){
  // returns from 0 to 1
  double temp = (s - t) / (h_n*sqrt(2.0));
  return (erf(temp) * 0.5 + 0.5 );
}

//[[Rcpp::export]]
double KerWgt(double s, double t, double h_n){
  double pi=3.14159265;
  return(1/sqrt(2.0*pi*pow(h_n,2))*exp(-pow(s-t,2)/(2*pow(h_n,2))));
}


//[[Rcpp::export]]
double M_ker(NumericVector M_seq, NumericVector t_seq, double s, double h_n){
  double WgtSum = 0;
  double ttlWgt = 0;
  double n = M_seq.size();
  for (int i=0; i<n; i++){
    double weight = KerWgt(t_seq[i], s, h_n);
    WgtSum += weight*M_seq[i];
    ttlWgt += weight;
  }
  return(WgtSum/ttlWgt);
  
}


// [[Rcpp::export]]
List data_ker_C_v3(NumericVector id, NumericVector Y, IntegerVector delta, 
                   NumericVector X1, NumericVector X2, NumericVector t,
                   List list_for_t_star, double h_n,
                   int n_star, NumericVector id_u){
  // asumme "list_for_t_star" is ordered by id
  
  List out;
  
  NumericVector id_star(n_star);
  IntegerVector delta_star(n_star);
  NumericVector X1_star(n_star);
  NumericVector X2_star(n_star);
  NumericVector t_star(n_star);
  
  int counter = 0;
  for (int i=0; i<list_for_t_star.size(); i++){
    NumericVector t_star_i = list_for_t_star[i];
    
    int id_i = id_u[i];
    
    NumericVector X1_i = X1[id==id_i];                             // the available X measurements
    NumericVector X2_i = X2[id==id_i];     
    NumericVector t_i = t[id==id_i];                               // the t at which Xs are measured
    
    
    for (int j=0; j<t_star_i.size(); j++){
      double t_j = t_star_i[j];
      double X1_ker = M_ker(X1_i, t_i, t_j, h_n);
      double X2_ker = M_ker(X2_i, t_i, t_j, h_n);
      
      X1_star[counter] = X1_ker;
      X2_star[counter] = X2_ker;
      t_star[counter] = t_j;
      id_star[counter] = id_i;
      IntegerVector temp = delta[id==id_i];
      delta_star[counter] = temp[0];
      counter ++;
    }
  }
  
  out["id_star"] = id_star;
  out["delta_star"] = delta_star;
  out["X1_star"] = X1_star;
  out["X2_star"] = X2_star;
  out["t_star"] = t_star;
  
  
  return(out);
}  


// [[Rcpp::export]]
List gen_index_C_v2(IntegerVector dat_id, NumericVector dat_delta, NumericVector dat_Y,
                    IntegerVector id, NumericVector t, NumericVector delta, NumericVector X1,
                    NumericVector X2, IntegerVector group){
  // group = 1 returned from data_ker_C_v2 means original data
  
  NumericVector d1 = delta[(delta==1) & (group==1)];
  List index(d1.size());
  
  IntegerVector id_ev = dat_id[dat_delta == 1];   //unique ids of subjects who had events
  
  int counter = 0;
  for (int i=0; i<id_ev.size(); i++){
    int id_i = id_ev[i];
    
    IntegerVector v = gen_seq_C(0, id.size()-1);
    IntegerVector idx_i = v[(id==id_i) & (group==1)];   // find index of rows of the subject in the combined data set
    
    NumericVector time_i = dat_Y[dat_id==id_i];                 // event time of subject i
    IntegerVector at_risk_id = dat_id[(dat_Y>=time_i[0]) & (dat_id!=id_i)];        // find the id of at risk population excluding i
    
    for (int j=0; j<idx_i.size(); j++){
      IntegerVector temp(at_risk_id.size()+1);
      int idx_ij = idx_i[j];
      temp[0] = idx_ij + 1;                           // store the index of jth measurement for subject i; +1 for R
      double measure_time = t[idx_ij];                // measurement time of the subject i
      
      int counter2 = 1;
      for (int k=0; k<at_risk_id.size(); k++){
        int id_k = at_risk_id[k];
        IntegerVector row = v[(id==id_k) & (t==measure_time)];
        temp[counter2] = row[0]+1;               // +1 for R            
        counter2++;
      }
      
      index[counter] = temp;
      counter++;
    }
    
  }
  
  return index;
  
}
// separte generate kernel smoothing data and generate index to make sure kernel smoothed data do nor enter the index of rows of subjects who had events;



// [[Rcpp::export]]
double pos_like_two_markers_ker_C(NumericVector coefs, int no_marker, int df, NumericVector X1, NumericVector X2, NumericVector t, List index, double h_n){
  
  double loglike = 0;
  
  for (int i = 0; i < index.size(); i++){
    
    NumericVector idx = index[i];
    NumericVector score(idx.size());                 // create a vector of the same length as idx
    for (int j = 0; j < idx.size(); j++){
      double t_ij = t[idx[0]-1];
      
      NumericVector t_ij_frac = fracpoly(t_ij);
      
      double beta1 = 0;
      for (int z = 0; z < df; z++){
        if (z==0) { 
          beta1 += t_ij_frac[z]*1;
        }else{
          beta1 += t_ij_frac[z]*coefs[z-1];
        }
      }
      
      double beta2 = 0;
      for (int z = 0; z < df; z++){
        beta2 += t_ij_frac[z]*coefs[z+df-1];
      }
      
      score[j] = beta1*X1[idx[j]-1] + beta2*X2[idx[j]-1];
    }
    
    double numer = 0;
    for (int k = 0; k < score.size(); k++){
      numer += KerSmooth(score[0], score[k], h_n);    // smoothed using gaussian kernel
    }
    double loglike_ij =  log(numer/score.size() );
    
    loglike += loglike_ij;
  }
  
  
  return -loglike;
}
// 20200409update: p-1 coefficients is correct; p is wrong
// 20200506update: no need to add a tiny number in log

// [[Rcpp::export]]
double neg_like_two_markers_ker_C(NumericVector coefs, int no_marker, int df, NumericVector X1, NumericVector X2, NumericVector t, List index, double h_n){

  
  double loglike = 0;
  
  for (int i = 0; i < index.size(); i++){
    
    NumericVector idx = index[i];
    NumericVector score(idx.size());                 // create a vector of the same length as idx
    for (int j = 0; j < idx.size(); j++){
      double t_ij = t[idx[0]-1];
      
      NumericVector t_ij_frac = fracpoly(t_ij);
      
      double beta1 = 0;
      for (int z = 0; z < df; z++){
        if (z==0) { 
          beta1 += t_ij_frac[z]*(-1);
        }else{
          beta1 += t_ij_frac[z]*coefs[z-1];
        }
      }
      
      double beta2 = 0;
      for (int z = 0; z < df; z++){
        beta2 += t_ij_frac[z]*coefs[z+df-1];
      }
      score[j] = beta1*X1[idx[j]-1] + beta2*X2[idx[j]-1];
    }
    
    double numer = 0;
    for (int k = 0; k < score.size(); k++){
      numer += KerSmooth(score[0], score[k], h_n);    // smoothed using gaussian kernel
    }
    double loglike_ij =  log(numer/score.size());
    
    loglike += loglike_ij;
  }
  
  
  return -loglike;
}
// 20200409 update: p-1 coefficients is correct. p is wrong.
// 20200506 update: no need to add a tiny number in log
