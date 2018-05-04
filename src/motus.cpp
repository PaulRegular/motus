// Calculates the negative log likliehood of a one behaviour DCRW.
// Authors: Marie Auger-Methe, Ian Jonsen, Christoffer M. Albertsen
/* Code associated with case study 1 from:
 * Template Model Builder: a promising tool for modelling the movement of marine animals.
 * Auger-Methe, M, Albertsen, CM, Jonsen, ID, Derocher, AE, Lidgard, D, Studholme, KR,
 * Crossin, GT, Bowen, WD, Mills Flemming J.
 */

// modified by Paul Regular to fit to accoustic telemetry data

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

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Input data
    DATA_VECTOR(y_lon);       // Observed locations
    DATA_VECTOR(y_lat);
    DATA_VECTOR(delta_t);     // Time interval between two consecutive location
    DATA_VECTOR(obs_sd_lon);  // normal distribution observation error parameters from obs_error model
    DATA_VECTOR(obs_sd_lat);
    DATA_VECTOR(tau_lon);     // t-distribution observation error parameters from ""
    DATA_VECTOR(tau_lat);
    DATA_VECTOR(nu_lon);
    DATA_VECTOR(nu_lat);
    DATA_VECTOR(scale_lon);   // Cauchy distribution observation error parameters from ""
    DATA_VECTOR(scale_lat);
    DATA_INTEGER(n);          // Number of records
    DATA_INTEGER(dist);       // Switch for the distribution to use for observation error

    // Input parameters
    PARAMETER_VECTOR(logit_gamma);   // Autocorrelation - logit because 0 < gamma < 1
    PARAMETER(log_sd_gamma);         // Error for gamma random walk
    PARAMETER(log_sd_lon);           // Process error (sd) in lon - log because sd > 0
    PARAMETER(log_sd_lat);           // Process error (sd) in lat
    PARAMETER(log_alpha_lon);        // Proportional constant for observation error > 0 (because sd > 0)
    PARAMETER(log_alpha_lat);

    // The unobserved locations of the animal, i.e states
    PARAMETER_VECTOR(x_lon_km);
    PARAMETER_VECTOR(x_lat_km);

    // Transformation of the input parameters to model format
    /* These transformations are made to insured that the parameters have sensical values.
     They do not change the model, they are only a computational trick. */
    vector<Type> gamma = 1.0 / (1.0 + exp(-logit_gamma));      // logit-1 b/c we want 0 < gamma < 1
    Type sd_gamma = exp(log_sd_gamma);
    Type sd_lon = exp(log_sd_lon);                             // exp-log b/c we want sd > 0
    Type sd_lat = exp(log_sd_lat);
    Type alpha_lon = exp(log_alpha_lon);
    Type alpha_lat = exp(log_alpha_lat);
    vector<Type> x_lon = x_lon_km * Type(1000.0);
    vector<Type> x_lat = x_lat_km * Type(1000.0);

    // Create a variable that will keep track of the negative log likelihood (nll)
    parallel_accumulator<Type> nll(this);

    // Create a temporary variables to be used below
    Type tmp_gamma;
    Type tmp_lon;
    Type tmp_lat;

    // Process equation
    for (int i = 1; i < n; ++i) {
        nll -= dnorm(logit_gamma(i), logit_gamma(i - 1), delta_t(i) * sd_gamma, true);   // Assume gamma follows a random walk
        if (i == 1) {
            nll -= dnorm(x_lon(i), x_lon(i - 1), delta_t(i) * sd_lon, true);             // Assume a simple random walk from the first location
            nll -= dnorm(x_lat(i), x_lat(i - 1), delta_t(i) * sd_lat, true);
        }
        if (i >= 2) {
            tmp_gamma = 1.0 / (1.0 + exp(-logit_gamma(i)));
            tmp_lon = x_lon(i - 1) + tmp_gamma * (delta_t(i)/delta_t(i - 1)) * (x_lon(i - 1) - x_lon(i - 2));
            tmp_lat = x_lat(i - 1) + tmp_gamma * (delta_t(i)/delta_t(i - 1)) * (x_lat(i - 1) - x_lat(i - 2));
            nll -= dnorm(x_lon(i), tmp_lon, delta_t(i) * sd_lon, true);
            nll -= dnorm(x_lat(i), tmp_lat, delta_t(i) * sd_lat, true);
        }
    }

    // Observation equation
    for (int i = 0; i < n; ++i) {
        if (dist == 0) {
            nll -= dnorm(y_lon(i), x_lon(i), alpha_lon * obs_sd_lon(i), true);
            nll -= dnorm(y_lat(i), x_lat(i), alpha_lat * obs_sd_lat(i), true);
        }
        if (dist == 1) {
            tmp_lon = (y_lon(i) - x_lon(i)) / (alpha_lon * tau_lon(i));
            tmp_lat = (y_lat(i) - x_lat(i)) / (alpha_lat * tau_lat(i));
            nll -=  log(1 / (alpha_lon * tau_lon(i))) + dt(tmp_lon, nu_lon(i), true);
            nll -=  log(1 / (alpha_lat * tau_lat(i))) + dt(tmp_lat, nu_lat(i), true);
        }
        if (dist == 2) {
            nll -= dcauchy(y_lon(i), x_lon(i), alpha_lon * scale_lon(i), true);
            nll -= dcauchy(y_lat(i), x_lat(i), alpha_lat * scale_lat(i), true);
        }
    }

    ADREPORT(x_lon);
    ADREPORT(x_lat);

    return nll;

}
