/* 
* Roboception GmbH 
* Munich, Germany 
* www.roboception.com 
* 
* Copyright (c) 2019 Roboception GmbH 
* All rights reserved 
* 
* Author: Monika Florek-Jasinska
*/
#ifndef PARTICLE_FILTER_MULTIVARIATEPROB_H
#define PARTICLE_FILTER_MULTIVARIATEPROB_H

namespace gauss
{
class MultivGauss
{
  public:
    MultivGauss(double sig_x, double sig_y) : sig_x_(sig_x), sig_y_(sig_y) {};

    double multiv_prob(double x_dist, double y_dist);


  private:
    double sig_x_, sig_y_;

};
}

#endif //PARTICLE_FILTER_MULTIVARIATEPROB_H
