#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Input data
  DATA_VECTOR(lon);
  DATA_VECTOR(lat);
  DATA_IVECTOR(n_i);
  
  // Parameters
  PARAMETER_VECTOR(log_tau_lon);
  PARAMETER_VECTOR(log_nu_lon);
  PARAMETER_VECTOR(log_tau_lat);
  PARAMETER_VECTOR(log_nu_lat);
  Type mu = 0.0;
  
  // Transformations
  vector<Type> tau_lon = exp(log_tau_lon);
  vector<Type> nu_lon = exp(log_nu_lon);
  vector<Type> tau_lat = exp(log_tau_lat);
  vector<Type> nu_lat = exp(log_nu_lat);
  
  // Model
  parallel_accumulator<Type> nll(this);
  for(int i = 0; i < lon.size(); ++i){
    nll -= dt((lon(i) - mu)/tau_lon(n_i(i)), nu_lon(n_i(i)), true) - log(tau_lon(n_i(i)));
    nll -= dt((lat(i) - mu)/tau_lat(n_i(i)), nu_lat(n_i(i)), true) - log(tau_lat(n_i(i)));
  }
  
  // Reports
  REPORT(tau_lon);
  REPORT(nu_lon);
  REPORT(tau_lat);
  REPORT(nu_lat);
  return nll;
  
}
