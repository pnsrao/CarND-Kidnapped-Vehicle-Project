/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

bool open_loop = false;
bool genie_pos = false;
double gt_x;
double gt_y;
double gt_theta;
int time_step = -1;

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  time_step++;
  gt_x = x;gt_y = y;gt_theta = theta;
  cout<<time_step<<" gt "<<gt_x<<" "<<gt_y<<" "<<gt_theta<<endl;
  if(open_loop || genie_pos){
    if(!is_initialized){
      num_particles = 1;
      Particle p;p.x = x;p.y = y;p.weight =  1;p.theta = theta;
      particles.clear();
      particles.push_back(p);
      is_initialized = true;
      return;
    }
    else if(genie_pos){
      particles[0].x = x;particles[0].y = y;particles[0].theta= theta;
      return;
    }
    else{
      return;
    }
  }
  if(is_initialized){return;}
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
        default_random_engine gen;
  	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);
	
	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	num_particles = 10;
	for(int ii=0;ii<num_particles;ii++){
	  Particle p;
	  p.x = dist_x(gen);
	  p.y = dist_y(gen);
	  p.theta = dist_theta(gen);
	  p.weight = 1;
	  particles.push_back(p);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  if(genie_pos){return;}
  bool verbose = false;
        default_random_engine gen;
  	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(0, std_pos[0]);
	
	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
        for(int ii=0;ii<num_particles;ii++){
	if(verbose){cout <<"Pre "<<particles[ii].x<<" "<<particles[ii].y<<endl;}
	  float temp_theta = particles[ii].theta;
	  if(yaw_rate){
	    particles[ii].x += (velocity/yaw_rate)*(sin(temp_theta+yaw_rate*delta_t)-sin(temp_theta));
	  }
	  else{
	    particles[ii].x += velocity*cos(temp_theta)*delta_t;
	  }
	  if(!open_loop){particles[ii].x += dist_x(gen);}
	  if(yaw_rate){
	    particles[ii].y += (velocity/yaw_rate)*(cos(temp_theta)-cos(temp_theta+yaw_rate*delta_t));
	  }
	  else{
	    particles[ii].y += velocity*sin(temp_theta)*delta_t;
	  }
	  if(!open_loop){particles[ii].y += dist_y(gen);}
	  particles[ii].theta += yaw_rate*delta_t;
	  if(!open_loop){particles[ii].theta += dist_theta(gen);}
	if(verbose){cout <<"Prediction done "<<particles[ii].x<<" "<<particles[ii].y<<endl;}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  bool verbose = false;  
  if(verbose){cout<<observations.size()<<" "<<predicted.size()<<endl;}
        for(int ii=0;ii<observations.size();ii++){
	  float min_distance = sqrt(pow((predicted[0].x - observations[ii].x),2)+pow((predicted[0].y - observations[ii].y),2));
	  int min_index = 0;
	  for(int jj=1;jj<predicted.size();jj++){
	    float mydistance = sqrt(pow((predicted[jj].x - observations[ii].x),2)+pow((predicted[jj].y - observations[ii].y),2));
	    if(mydistance < min_distance){
	      min_distance = mydistance;
	      min_index = jj;
	    }
	  }
	  observations[ii].id = predicted[min_index].id;
	  // Print distances for debug
	  if(false){
	    cout <<observations[ii].x<<" "<<observations[ii].y<<" "<<predicted[min_index].x<<" "<<predicted[min_index].y<<endl;
	  }
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  bool verbose = false;if(verbose){cout<<"In updateWeights "<<num_particles<<endl;}
        for(int ii=0;ii<num_particles;ii++){
	  vector<LandmarkObs> this_obs(observations.size());
	  Particle p = particles[ii];
	  for(int jj=0;jj<observations.size();jj++){
	    //Transformations
	    this_obs[jj].x = p.x + cos(p.theta)*observations[jj].x - sin(p.theta)*observations[jj].y;
	    this_obs[jj].y = p.y + sin(p.theta)*observations[jj].x + cos(p.theta)*observations[jj].y;
	  }
	  if(verbose){cout<<"Transformation done"<<endl;}
	  int num_landmarks = map_landmarks.landmark_list.size();
	  vector<LandmarkObs> preds;
	  for(int jj=0;jj<num_landmarks;jj++){
	    //Populate structure for nearest neighbor association
	    // Check sensor range
	    float lx = map_landmarks.landmark_list[jj].x_f;
	    float ly = map_landmarks.landmark_list[jj].y_f;
	    if((fabs(lx-p.x) <= sensor_range) && (fabs(ly-p.y)<=sensor_range)){
	      LandmarkObs temp;
	      // Index landmarks by the index in the vector, not by the id in the file.
	      temp.id = jj; //map_landmarks.landmark_list[jj].id_i;
	      temp.x = lx;
	      temp.y = ly;
	      preds.push_back(temp);
	    }
	  }
	  if(preds.size() == 0){
	    particles[ii].weight = 1e-15;
	    this_obs.clear();
	    return;
	  }
	  if(verbose){cout<<"Association prep done"<<endl;}
	  dataAssociation(preds,this_obs);
	  if(verbose){cout<<"Association done"<<endl;}
	  float prob_meas = 1;
	  for(int jj=0;jj<observations.size();jj++){
	    //Find associated landmark id
	    int l_id = this_obs[jj].id;
	    float myexp = pow((this_obs[jj].x - map_landmarks.landmark_list[l_id].x_f)/std_landmark[0],2)+pow((this_obs[jj].y - map_landmarks.landmark_list[l_id].y_f)/std_landmark[1],2);
	    float this_prob = exp(-myexp/2.0)/(2*M_PI*std_landmark[0]*std_landmark[1]);
	    prob_meas *= this_prob;
	    if(verbose){
	      //cout<<ii<<" "<<jj<<" "<<this_obs[jj].x<<" "<<map_landmarks.landmark_list[l_id].x_f<<" "<<myexp<<endl;
	      cout<<ii<<" "<<jj<<" "<<prob_meas<<" "<<this_prob<<" "<<myexp<<endl;
	    }
	  }
	  particles[ii].weight = prob_meas;
	  this_obs.clear();
	  preds.clear();
	  if(verbose){cout<<"Weights done"<<endl;}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  if(genie_pos || open_loop){return;}
  vector<double> myweights(num_particles);
  default_random_engine gen;
  for(int ii=0;ii<num_particles;ii++){
    myweights[ii] = particles[ii].weight;
  }
  discrete_distribution<> dist_new(myweights.begin(),myweights.end());
  vector<Particle> temp_particles = particles;
  for(int ii=0;ii<num_particles;ii++){
    int index = dist_new(gen);
    particles[ii] = temp_particles[index];
  }
  temp_particles.clear();
  

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
