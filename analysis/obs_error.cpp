#include <TMB.hpp>

template<class Type>
Type dcauchy(Type x, Type location, Type scale, int give_log){
    Type pi = 3.1415926535897932384;
    Type y = (x - location) / scale;
    Type res = Type(1.0) / (pi * scale * (Type(1.0) + y * y));
    if(give_log) {
        return(log(res));
    } else {
        return(res);
    }
}
VECTORIZE4_ttti(dcauchy);

// Modification of Anders Nielsen's code
template<class Type>
Type drobust(Type x, Type mean, Type sd, Type p, Type df, int give_log = 0){
    Type z = (x - mean) / sd;
    Type logres = log(Type(1.0) / sd * ((1.0 - p) *
        dnorm(z, Type(0.0), Type(1.0), false) + p * dt(z, df, false)));
    if(!give_log){
        return exp(logres);
    }else{
        return logres;
    }
}
VECTORIZE6_ttttti(drobust);


template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(lon);
  DATA_VECTOR(lat);
  DATA_IVECTOR(n_i);
  DATA_INTEGER(dist); // toggle for distribution to fit

  Type mu = 0.0;

  parallel_accumulator<Type> nll(this);

  if (dist == 0) { // normal distribution

      PARAMETER_VECTOR(log_sd_lon);
      PARAMETER_VECTOR(log_sd_lat);

      vector<Type> sd_lon = exp(log_sd_lon);
      vector<Type> sd_lat = exp(log_sd_lat);

      for (int i = 0; i < lon.size(); ++i){
          nll -= dnorm(lon(i), mu, sd_lon(n_i(i)), true);
          nll -= dnorm(lat(i), mu, sd_lat(n_i(i)), true);
      }

      REPORT(sd_lon);
      REPORT(sd_lat);

  }

  if (dist == 1) { // t-distribution

      PARAMETER_VECTOR(log_tau_lon);
      PARAMETER_VECTOR(log_nu_lon);
      PARAMETER_VECTOR(log_tau_lat);
      PARAMETER_VECTOR(log_nu_lat);

      vector<Type> tau_lon = exp(log_tau_lon);
      vector<Type> nu_lon = exp(log_nu_lon);
      vector<Type> tau_lat = exp(log_tau_lat);
      vector<Type> nu_lat = exp(log_nu_lat);

      for (int i = 0; i < lon.size(); ++i){
          nll -= dt((lon(i) - mu)/tau_lon(n_i(i)), nu_lon(n_i(i)), true) - log(tau_lon(n_i(i)));
          nll -= dt((lat(i) - mu)/tau_lat(n_i(i)), nu_lat(n_i(i)), true) - log(tau_lat(n_i(i)));
      }

      REPORT(tau_lon);
      REPORT(nu_lon);
      REPORT(tau_lat);
      REPORT(nu_lat);

  }

  if (dist == 2) { // Cauchy distribution

      PARAMETER_VECTOR(log_scale_lon);
      PARAMETER_VECTOR(log_scale_lat);

      vector<Type> scale_lon = exp(log_scale_lon);
      vector<Type> scale_lat = exp(log_scale_lat);

      for(int i = 0; i < lon.size(); ++i){
          nll -= dcauchy(lon(i), mu, scale_lon(n_i(i)), true);
          nll -= dcauchy(lat(i), mu, scale_lat(n_i(i)), true);
      }

      REPORT(scale_lon);
      REPORT(scale_lat);

  }

  if (dist == 3) { // normal - t mixture

      PARAMETER_VECTOR(log_sd_lon);
      PARAMETER_VECTOR(log_sd_lat);
      PARAMETER_VECTOR(logit_p); // proportion t
      PARAMETER_VECTOR(log_df);  // degrees of freedom

      vector<Type> sd_lon = exp(log_sd_lon);
      vector<Type> sd_lat = exp(log_sd_lat);
      vector<Type> p = 1.0 / (1.0 + exp(-logit_p));
      vector<Type> df = exp(log_df);

      for(int i = 0; i < lon.size(); ++i){
          nll -= drobust(lon(i), mu, sd_lon(n_i(i)), p(n_i(i)), df(n_i(i)), true);
          nll -= drobust(lat(i), mu, sd_lat(n_i(i)), p(n_i(i)), df(n_i(i)), true);
      }

      REPORT(sd_lon);
      REPORT(sd_lat);
      REPORT(p);
      REPORT(df);

  }

  return nll;

}
