/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include <tuple>
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  particles.resize(num_particles);
  std::random_device rd{};
  std::mt19937 gen{rd()};

  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> x_norm{x, std[0]};
  std::normal_distribution<> y_norm{y, std[1]};
  std::normal_distribution<> theta_norm{theta, std[2]};

  for (auto &particle : particles)
  {

  particle.x = x + x_norm(gen);
  particle.y = y + y_norm(gen);
  particle.theta = theta + theta_norm(gen);
    particle.weight = 1;
    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::random_device rd{};
  std::mt19937 gen{rd()};
delta_t = delta_t/100000;
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean

  for (auto &particle : particles)
  {
    std::normal_distribution<> x{particle.x, std_pos[0]};
    std::normal_distribution<> y{particle.y, std_pos[1]};
    std::normal_distribution<> theta{particle.theta, std_pos[2]};

    particle.x = particle.x + velocity * sin(yaw_rate) * delta_t;// + x(gen); //TODO possible wrong
    particle.y = particle.y + velocity * cos(yaw_rate) * delta_t;// + y(gen);
    particle.theta = particle.theta + theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  auto convertToMap = [](vector<LandmarkObs> observations, const Particle &particle)
  {
    for (auto &obs: observations)
    {
      //double vec = sqrt(pow(observation.x, 2) + pow(observation.y , 2));
//      observation.x = particle.x + vec*sin(particle.theta); //TODO if correct?
//      observation.y = particle.y + vec*cos(particle.theta); //TODO if correct?
      double cos_theta = cos(particle.theta);
      double sin_theta = sin(particle.theta);
      double obs_x = particle.x + obs.x * cos_theta - obs.y * sin_theta;
      double obs_y = particle.y + obs.x * sin_theta + obs.y * cos_theta;
      obs.x = obs_x;
      obs.y = obs_y;

    }
    return observations;
  };

//  auto getPredicted = [&](const Particle &particle)
//  {
//    vector<LandmarkObs> map_obs;
//    for (auto &obs: map_landmarks.landmark_list)
//    {
//      LandmarkObs observation;
//      double vec = sqrt(pow(obs.x_f - particle.x, 2) + pow(obs.y_f - particle.y, 2));
//      observation.x = vec * sin(particle.theta); //TODO if correct?
//      observation.y = vec * cos(particle.theta); //TODO if correct?
//      map_obs.push_back(observation);
//    }
//    return map_obs;
//  };
  auto getClosestDistXY = [&](const LandmarkObs &obs_in_map)
  {
//    auto distance_xy = std::make_tuple(sensor_range, sensor_range);
//    double closest = sensor_range;
//    for (auto obs: observations)
//    {
//      double dist_x =pow(obs.x - map_landmark.x, 2);
//      double dist_y = pow(obs.y - map_landmark.y, 2);
//      double res = sqrt( dist_x + dist_y);
//      if (res < closest)
//      {
//        closest = res;
//        distance_xy = std::make_tuple(dist_x, dist_y);
//      }
//
//    }
//    return distance_xy;
    auto distance_xy = std::make_tuple(sensor_range, sensor_range);
    double closest = pow(sensor_range, 2);
    for (auto obs : map_landmarks.landmark_list)
    {
      double dist_x =pow(obs.x_f - obs_in_map.x, 2);
      double dist_y = pow(obs.y_f - obs_in_map.y, 2);
      double res = sqrt(dist_x + dist_y);
      if (res < closest)
      {
        closest = res;
        distance_xy = std::make_tuple(dist_x, dist_y);
      }

    }
    return distance_xy;
  };
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double sum_of_weights = 0;
  gauss::MultivGauss multiv_gauss(std_landmark[0], std_landmark[1]);
  for (auto &particle : particles)
  {
    //auto map_landmarks = getPredicted(particle);
    particle.weight = 1;
    auto observations_in_map = convertToMap(observations, particle);
    for (auto obs : observations_in_map){
      auto dist = getClosestDistXY(obs);
      double prob = multiv_gauss.multiv_prob(std::get<0>(dist),std::get<1>(dist));
      if (prob < 0.0001) prob =0.0001;
      particle.weight = particle.weight * prob;
    }
    std::cout<<particle.weight<<std::endl;

  //sum_of_weights += particle.weight;
  }

  for (auto &particle : particles) particle.weight = particle.weight / sum_of_weights;
}

void ParticleFilter::resample()
{
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::default_random_engine generator;
  std::vector<double> weights(num_particles);
  for (int i = 0; i < num_particles; i++){
    weights[i] = particles[i].weight;
  }

  std::discrete_distribution<int> distribution(weights.begin(), weights.end());

  std::vector<Particle> new_particles(num_particles);

  for (int i=0; i<num_particles; ++i) {
    int number = distribution(generator);
    new_particles[i] = particles[number];
  }


}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
