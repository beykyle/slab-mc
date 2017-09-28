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
    Particle(float x , float mu);
    ~Particle(void);
};

Particle::Particle()
{
  x = 0;
  mu = 0;
  alive = true;
};

Particle::Particle(float x_in , float mu_in) 
{
  x = x_in;
  mu = mu_in;
  alive = true;
};

class Slab {
  public:
    double               width;
    int                  nCells;
    std::vector <double> scatter_xs;
    std::vector <double> total_xs;
    std::vector <double> fission_xs;
};


class Leakage {
  public:
    bool     vector;
    bool     total;
    int      nbins;

    void     incrementCount(void);
    void     setVal(double nHistories , Slab s);
  private:
    double   count;
    double   squared;
    double   hist;
    double   value;
    double   uncertainty;
};

void Leakage::setVal(double nHistories , Slab s) { 
  // Finalizes the leakage tally and calculates uncertainty
  value =  count / nHistories;
  uncertainty = sqrt( ( squared / nHistories - pow(value,2) ) / nHistories );
};

void Leakage::incrementCount(void) {
  count += 1.0;
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
    double xi    = Urand();
    bool   found = false;
    float  prob = 0;
    float  mu;
    
    std::cout << "xi is " << xi << std::endl;
    for(int i = 0; i < s.mu_probs.size(); ++i) {
        if(xi > prob) {
          std::cout << "prob is " << prob << std::endl;
          prob += s.mu_probs[i];
        }
        else if(found == false) {
          mu = s.mu_bins[i];
          found = true;
        }
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

void transport(Tally &leak_tally , Tally &flux_tally, Slab s) {
//  this function takes in pointers to a slab and tally objects, and a pointer to a stack of particles, 
//  modifies tally objects tally objects 
//  for each particle in the stack
//  run transport
    // Particle * tmp = stack.front();
    // stack.pop() // get rid of the pointer to the particle 
    // delete(&tmp) // destruct the memory holding the particle
};

// ********************************************************************************************************* //
//  Main Function: creates input and tally objects, calls and times transport(), calls output functions      //
//												             //
// ********************************************************************************************************* //

int main()  {
  
  // set up transport problem
  //TODO create input class constructor that parses from input file
  Input I; 
  I.nHistories       = 1;
  I.slab.width       = 4.0;
  I.slab.nCells      = 4;
  I.slab.total_xs    = {0.0 , 0.3 , 1.0 , 1.0 } ;
  I.slab.scatter_xs  = {0.0 , 0.25 , 0.7 , 0.1};
  I.slab.fission_xs  = {0.0, 0.0 , 0.0 , 0.0};

  // define a source
  Source s;
  s.location =   2.0;
  s.mu_bins  = {-1.0 , -0.8 , -0.6, -0.4 , 0.0 , 0.8};
  s.mu_probs = { 0.1 ,  0.2 , 0.2 , 0.1 ,  0.2 , 0.2};

  // initialize a particle stack
  std::queue<Particle*> stack; 
  stack = makeStack(s , I.nHistories);
  std::cout << " the angle of the first particle on the stack isssssss: "<< stack.front()->mu << std::endl;
  // create and populate a dictionary of Tally objects
  // with the key being the quantity being tallied
  Leakage leak_tally;
  Flux    flux_tally;
  flux_tally.nbins = 1000; //the mesh for the flux tally has 1000 bins

  //I.tallies["Leakage Probability"] = leak_tally;
  //I.tallies["Flux"] = flux_tally;
  

  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  // transport(I.tallies["Leakage Probability"], I.tallies["Flux"], I.slab);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  
  // output results 
  // terminal_out( I.tallies , duration);
  // file_out(I.slab , I.tallies , duration); 
  
  return(0);
}
