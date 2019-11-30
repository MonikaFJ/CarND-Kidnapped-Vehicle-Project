#include "multiv_gauss.h"
#include <cmath>

namespace gauss
{
double MultivGauss::multiv_prob(double x_obs, double y_obs)
{
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x_ * sig_y_);

  // calculate exponent
  double exponent;
  exponent = x_obs / (2 * pow(sig_x_, 2))
             + y_obs / (2 * pow(sig_y_, 2));

  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);

  return weight;
}

}