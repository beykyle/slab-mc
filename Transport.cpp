#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <numeric>
#include <queue>
#include "Random.h"

//TODO make input class  in Input.cpp and Input.h, with read from file function
//TODO replace boost map with templated classes

class Source {
    public:
      double location;
      std::vector<double> mu_bins;  //discrete values for source emission angle
      std::vector<double> mu_probs; //probability for each value
};

class Particle {
  public:
    double x, mu;
    int    cell;
    bool   alive;
    Particle();
    Particle(double x , double mu);
    ~Particle(void);
};

Particle::Particle()
{
  x = 0;
  mu = 0;
  alive = true;
};

Particle::Particle(double x_in , double mu_in) 
{
  x = x_in;
  mu = mu_in;
  alive = true;
};

class Slab {
  public:
    double               width;
    int                  nCells;
    double scatter_xs;
    double total_xs;
    double fission_xs;
//    std::vector <double> scatter_xs;
//    std::vector <double> total_xs;
//    std::vector <double> fission_xs;
};


class Leakage {
  public:
    void       incrementCount(void);
    void       setVal(double nHistories);
    double     getVal(void);
               Leakage();
  private:
    double   count;
    double   squared;
    double   hist;
    double   value;
    double   uncertainty;
};

void Leakage::setVal(double nHistories ) { 
  // Finalizes the leakage tally and calculates uncertainty
  value =  count / nHistories;
  double squared = 1;
  uncertainty = sqrt( ( squared / nHistories - pow(value,2) ) / nHistories );
};

void Leakage::incrementCount(void) {
  count += 1.0;
};

Leakage::Leakage() {
  count = 0;
  squared = 0;
  hist = 0;
  value = 0;
  uncertainty = 0;
};

double Leakage::getVal(void) {
  return(value);
};


class Flux {
  public:
    bool                 vector;
    bool                 total;
    int                  nbins;

    void                 setVal(double nhistories , Slab s);
    void                 addLength(int i , double dist);
  private:
    double               value;
    double               uncertainty;
    std::vector <double> path_lengths;
    std::vector <double> values;
    std::vector <double> uncertainties;
};

void Flux::addLength(int i , double dist) {
  path_lengths[i] += dist;
};

void Flux::setVal(double nHistories , Slab s) {  
  // Finalizes the flux tally and calculates uncertainty
  double  sum_path_lengths = 0;

  for(int i = 0; i < path_lengths.size(); ++i) {
    values.push_back(path_lengths[i] / nHistories / (s.width / nbins ) );
    uncertainties.push_back(0.0); 
  }

  value = sum_path_lengths;
  uncertainty = 0.0; 
};

class Tally {
  public:
    bool                 vector;
    bool                 total;
    double               value , uncertainty;
    std::vector <double> values, uncertainties;
};

class Input {
  public:  
    Slab                           slab;
    int                            tallyBins;
    double                         nHistories;
    std::map<std::string , Tally>  tallies;

};

// ********************************************************************************************************* //
//  makeStack function: creates the initial particle stack according to the source definition                //
//	                   							         	                                                 //
// ********************************************************************************************************* //


std::queue<Particle*> makeStack(Source s , int size) {
  // takes in a pointer to a Source object and the number of histories
  // outputs a particle stack
  
  std::queue<Particle*> stack; // fite me irl
  for(int i=0 ; i < size ; ++i) {

    //  sample from emission angle
    double  xi    = Urand();
    bool    found = false;
    double  prob  = 0;
    double  mu;
    
    for(int i = 0; i < s.mu_probs.size(); ++i) {
        if(xi > prob) {
          prob += s.mu_probs[i];
        }
        else if(found == false) {
          mu = s.mu_bins[i];
          found = true;
        }
    }

    if(found == false) {
      mu = s.mu_bins[s.mu_bins.size()-1];
    }
    
    Particle *p = new Particle(s.location , mu);
    // push a pointer to p to the stack
    stack.push(p);
  }
  return(stack);
};


// ********************************************************************************************************* //
//  Transport function: Sets up transport problem according to Input object, modifies Tally objects          //
//	                   							         	                                                 //
// ********************************************************************************************************* //

void transport(Leakage &leak_tally , Flux &flux_tally, Slab s , std::queue<Particle*> stack) 
{
//  this function takes in pointers to a slab and tally objects, and a pointer to a stack of particles, 
//  modifies tally objects tally objects 
//  for each particle in the stack
//  run transport
  double b = s.width / flux_tally.nbins;

  while( ! stack.empty() ){  

      while ( stack.front()->alive ) {
        // do transport
        double optical_dist = -log( Urand() ); // optical distance to collison
        double x1           = stack.front()->x;

        stack.front()->x += optical_dist * stack.front()->mu / s.total_xs;
        double delta_x = stack.front()->x - x1;
 
        // check for leakage
        if ( stack.front()->x < 0.0 || stack.front()->x > s.width ) {
          // leaked out
          stack.front()->alive = false;
          if ( stack.front()->x < 0.0 ) {
            // update distance traveled and position to start of the slab
            stack.front()-> x = 0;
            optical_dist = s.total_xs * x1 / std::abs(stack.front()->mu); 
          // tally if leaked out of right side
          }
          if ( stack.front()->x > s.width ) {
            //update leakage tally
            leak_tally.incrementCount();
            // update distance traveled and position to end of the slab 
            stack.front()->x = s.width;
            optical_dist = s.total_xs * ( stack.front()->x - x1 ) / std::abs(stack.front()->mu);
          }
        }

        // update path flux tally in each spatial bin
        // indexing from 0
        int    oldBin  = floor(x1 / b);
        int    newBin  = floor(stack.front()->x / b);          
      
        // if the particle stayed in the same bin,
        // add the whole path length to the sum in that bin
        if ( oldBin == newBin ) {
        //   path_lengths[oldBin] += optical_dist / s.total_xs;
        }
        //  otherwise iterate through the bins the particle passed through, 
        //  incrementing them by the path length traversed in that bin
        else {
          //sort the bins and positions the particles moved between
          bool forward = newBin > oldBin; 
        
          int    smallBin   = (forward) ? (oldBin):(newBin);
          int    largeBin   = (forward) ? (newBin):(oldBin);
          double small_x    = (forward) ? (x1):(stack.front()->x);
          double large_x    = (forward) ? (stack.front()->x):(x1);
        }
        //    path_lengths[smallBin] += dist_to_collision * (( smallBin +1 ) * b - small_x ) / delta_x;

        // determine next transport step
        if (stack.front()->alive == true) {
          // have a collision
          if ( Urand() < s.scatter_xs / s.total_xs ) {
            // particle scatters
            stack.front()->mu = 2.0 * Urand() - 1.0;
          }
          else {
            // absorbed :(
            stack.front()->alive = false;
          }
      }
    } // end of current particle

    // destruct particle object and remove ptr from stack
    Particle * tmp = stack.front();
    stack.pop(); // get rid of the pointer to the particle 
    delete(&tmp); // destruct the memory holding the particle
  } // end of particle stack

};

// ********************************************************************************************************* //
//  Main Function: creates input and tally objects, calls and times transport(), calls output functions      //
//												             //
// ********************************************************************************************************* //

int main()  {
  
  // set up transport problem
  //TODO create input class constructor that parses from input file
  Input I; 
  I.nHistories       = 1e8;
  I.slab.width       = 4.0;
  I.slab.nCells      = 4;
  I.slab.total_xs    = 1.0;
  I.slab.scatter_xs  = 0.0;
  I.slab.fission_xs  = 0.0;

//  I.slab.total_xs    = {0.0 , 0.3  , 1.0 , 1.0} ;
//  I.slab.scatter_xs  = {0.0 , 0.25 , 0.7 , 0.1};
//  I.slab.fission_xs  = {0.0 , 0.0  , 0.0 , 0.0};

  // define a source
  Source s;
//  s.location =   2.0;
//  s.mu_bins  = {-1.0 , -0.8 , -0.6, -0.4 , 0.0 , 0.8};
//  s.mu_probs = { 0.1 ,  0.2 , 0.2 , 0.1 ,  0.2 , 0.2};
  s.location =   0.0;
  s.mu_bins  = {1.0 };
  s.mu_probs = { 1.0};

  // initialize a particle stack
  std::queue<Particle*> stack; 
  stack = makeStack(s , I.nHistories);

  // create and populate a dictionary of Tally objects
  // with the key being the quantity being tallied
  Leakage leak_tally;
  Flux    flux_tally;
  flux_tally.nbins = 1000; //the mesh for the flux tally has 1000 bins

  
  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  transport(leak_tally, flux_tally, I.slab , stack);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  
  
  leak_tally.setVal(I.nHistories);

  std::cout<< "leakage prob: " << leak_tally.getVal() << std::endl;
  // output results 
  // terminal_out( I.tallies , duration);
  // file_out(I.slab , I.tallies , duration); 
  
  return(0);
}
