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
namespace pf
{
void ParticleFilter::init(double x, double y, double theta, double std[])
{
  num_particles = 10;
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

  std::random_device rd{};
  std::mt19937 gen{rd()};

  double vel_dt = velocity * delta_t;
  double yaw_rate_dt = yaw_rate * delta_t;
  double vel_over_yaw = velocity / yaw_rate;
  for (auto &particle : particles)
  {
    if (yaw_rate == 0)
    {
      particle.x += cos(particle.theta) * vel_dt;
      particle.y += sin(particle.theta) * vel_dt;
    }
    else
    {
      particle.x += vel_over_yaw * (sin(particle.theta + yaw_rate_dt) - sin(particle.theta));
      particle.y += vel_over_yaw * (cos(particle.theta) - cos(particle.theta + yaw_rate_dt));
      particle.theta += yaw_rate_dt;
    }
    std::normal_distribution<> x{particle.x, std_pos[0]};
    std::normal_distribution<> y{particle.y, std_pos[1]};
    std::normal_distribution<> theta{particle.theta, std_pos[2]};
    particle.x = x(gen);
    particle.y = y(gen);
    particle.theta = theta(gen);
  }

}

std::vector<LandmarkObs>
ParticleFilter::getLandmarksInRange(double sensor_range, const Particle &particle, const Map &map_landmarks)
{
  std::vector<LandmarkObs> landmarks_in_range;
  for (auto landmark : map_landmarks.landmark_list)
  {
    if ((pow(landmark.x_f - particle.x, 2) + pow(landmark.y_f - particle.y, 2)) <
        pow(sensor_range, 2))
    {
      landmarks_in_range.push_back(LandmarkObs(landmark.id_i, landmark.x_f, landmark.y_f));
    }
  }

  return landmarks_in_range;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  auto convertToMap = [](vector<LandmarkObs> observations, const Particle &particle)
  {
    for (auto &obs: observations)
    {
      double cos_theta = cos(particle.theta);
      double sin_theta = sin(particle.theta);
      double obs_x = particle.x + obs.x * cos_theta - obs.y * sin_theta;
      double obs_y = particle.y + obs.x * sin_theta + obs.y * cos_theta;
      obs.x = obs_x;
      obs.y = obs_y;
    }
    return observations;
  };

  gauss::MultivGauss multiv_gauss(std_landmark[0], std_landmark[1]);
  for (int i = 0; i < num_particles; i++)
  {
    particles[i].weight = 1;
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
    particles[i].associations.clear();

    auto observations_in_map = convertToMap(observations, particles[i]);
    auto landmarks_in_range = getLandmarksInRange(sensor_range, particles[i], map_landmarks);
    for (auto obs : observations_in_map)
    {
      double dist_x_min, dist_y_min, closest = pow(sensor_range, 2);
      for (int j = 0; j < landmarks_in_range.size(); j++)
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

      if (closest < pow(sensor_range, 2)) //not a single match found within sensor_range distance
      {
        double prob = multiv_gauss.multiv_prob(dist_x_min, dist_y_min);
        particles[i].weight *= prob;
        particles[i].associations.push_back(obs.id);
        particles[i].sense_x.push_back(obs.x);
        particles[i].sense_y.push_back(obs.y);
      }
      else
      {
        particles[i].weight *= 0; //TODO potential pitfall with noisy data
      }
    }

    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample()
{
  std::discrete_distribution<int> resample_dist(weights.begin(), weights.end());
  std::default_random_engine gen;

  vector<Particle> resampled_particles;
  resampled_particles.resize(num_particles);
  for (auto &single_particle :resampled_particles)
  {
    single_particle = particles[resample_dist(gen)];
  }
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
}