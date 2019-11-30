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
  num_particles = 5;  // TODO: Set the number of particles
  particles.resize(num_particles);
  weights.resize(num_particles, 1);

  std::random_device rd{};
  std::mt19937 gen{rd()};

  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> x_norm{x, std[0]};
  std::normal_distribution<> y_norm{y, std[1]};
  std::normal_distribution<> theta_norm{theta, std[2]};

  for (auto &particle : particles)
  {

    particle.x = x_norm(gen);
    particle.y = y_norm(gen);
    particle.theta = theta_norm(gen);
    particle.weight = 1;
  }
  is_initialized = true;

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
  //delta_t= delta_t/10000;
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean

  for (auto &particle : particles)
  {
    if (yaw_rate == 0)
    {
      particle.x += velocity * cos(particle.theta) * delta_t;
      particle.y += velocity * sin(particle.theta) * delta_t;
    }
    else
    {
      particle.x += velocity / yaw_rate * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta));
      particle.y += velocity / yaw_rate * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
      particle.theta += yaw_rate * delta_t;
    }
    std::normal_distribution<> x{particle.x, std_pos[0]};
    std::normal_distribution<> y{particle.y, std_pos[1]};
    std::normal_distribution<> theta{particle.theta, std_pos[2]};
    particle.x = x(gen);
    particle.y = y(gen);
    particle.theta = theta(gen);
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
      //TODO tested and ok.
      double cos_theta = cos(particle.theta);
      double sin_theta = sin(particle.theta);
      double obs_x = particle.x + obs.x * cos_theta - obs.y * sin_theta;
      double obs_y = particle.y + obs.x * sin_theta + obs.y * cos_theta;
      obs.x = obs_x;
      obs.y = obs_y;

    }
    return observations;
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

  gauss::MultivGauss multiv_gauss(std_landmark[0], std_landmark[1]);
  for (int i = 0; i < num_particles; i++)
  {
    //auto map_landmarks = getPredicted(particle);
    particles[i].weight = 1;
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
    particles[i].associations.clear();
//    particles[i].associations.resize(observations.size());
//    particles[i].sense_x.resize(observations.size());
//    particles[i].sense_y.resize(observations.size());
    auto observations_in_map = convertToMap(observations, particles[i]);
    std::vector<LandmarkObs> landmarks_in_range;
    for (auto landmark : map_landmarks.landmark_list)
    {
      if ((pow(landmark.x_f - particles[i].x, 2) + pow(landmark.y_f - particles[i].y, 2)) <
          pow(sensor_range, 2))
      {

        LandmarkObs new_landmark;

        new_landmark.x = landmark.x_f;
        new_landmark.y = landmark.y_f;
        new_landmark.id = landmark.id_i;
        landmarks_in_range.push_back(new_landmark);
      }
    }

    for (auto obs : observations_in_map)
    {

      double dist_x_min, dist_y_min = 50;
      double closest = 2500;
      for (int j = 0; j<landmarks_in_range.size(); j++)
      {
        double dist_x = pow(obs.x - landmarks_in_range[j].x, 2);
        double dist_y = pow(obs.y - landmarks_in_range[j].y, 2);
        double res = dist_x + dist_y;
        if (res < closest)
        {
          closest = res;
          dist_x_min = dist_x;
          dist_y_min = dist_y;
          obs.id = landmarks_in_range[j].id;
        }
      }
      double prob = multiv_gauss.multiv_prob(dist_x_min, dist_y_min);
      //if (prob < 0.000000000000001) prob = 0.000000000000001;
      particles[i].weight = particles[i].weight * prob;
      if (closest<2500)
      {
        particles[i].associations.push_back(obs.id);
        particles[i].sense_x.push_back(obs.x);
        particles[i].sense_y.push_back(obs.y);
      }
    }


    //if (particles[i].weight == 1) particles[i].weight = 0; //no valid connections found
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample()
{
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::discrete_distribution<int> resample_dist(weights.begin(), weights.end());
  std::default_random_engine gen;

  //Construct a vector to hold resampled particles
  vector<Particle> resampled_particles;
  resampled_particles.resize(num_particles);
  for(int i = 0; i<num_particles; ++i){
    int number = resample_dist(gen);
    resampled_particles[i] = particles[number];
  }

  //Assign resampled particles to current particles list
  particles = resampled_particles;

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
