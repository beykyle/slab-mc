#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <numeric>
#include <stack>
#include "Random.h"

//TODO make input class  in Input.cpp and Input.h, with read from file function

// ********************************************************************************************************* //
//  Source class                                                                                             //
// ********************************************************************************************************* //

class Source {
    public:
      double location;
      std::vector<double> mu_bins;  //discrete values for source emission angle
      std::vector<double> mu_probs; //probability for each value
};

// ********************************************************************************************************* //
//  Particle class                                                                                           //
// ********************************************************************************************************* //

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

// ********************************************************************************************************* //
//  Slab class                                                                                               //
// ********************************************************************************************************* //

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

// ********************************************************************************************************* //
//  Leakage class, for tallies of leakage probability                                                        //
// ********************************************************************************************************* //

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

// ********************************************************************************************************* //
//  Flux class, for tallies of path length density                                                           //
// ********************************************************************************************************* //


class Flux {
  public:
    std::vector <double> path_lengths;
                         
                         Flux(int meshBins , double x0 , double x1);
    std::vector<double>  getBounds(void);
    int                  getBins(void);
    double               getWidth(void);
    std::vector<double>  getVals(void);
    void                 setVals(double nhistories);
  private:
    double               width;
    int                  nbins;
    std::vector<double>  bounds;
    std::vector <double> values;
    std::vector <double> uncertainties;
};

int Flux::getBins(void) {
  return(nbins);
};

std::vector<double> Flux::getBounds(void) {
  return(bounds);
};

std::vector<double> Flux::getVals(void) {
  return(values);
}

Flux::Flux(int meshBins , double x0 , double x1) {
    nbins =  meshBins;
    path_lengths.resize(meshBins);
    width = x1 - x0;
    bounds = {x0 , x1};
};

void Flux::setVals(double nHistories ) {  
  // Finalizes the flux tally and calculates uncertainty
  for(int i = 0; i < path_lengths.size(); ++i) {
    values.push_back(path_lengths[i] / nHistories / (width / nbins ) );
    uncertainties.push_back(0.0); 
  }
};

// ********************************************************************************************************* //
//  Input class                                                                                              //
// ********************************************************************************************************* //

class Input {
  public:  
    Slab                           slab;
    double                         nHistories;
};

// ********************************************************************************************************* //
//  makeStack function: creates the initial particle stack according to the source definition                //
//	                   							         	                                                 //
// ********************************************************************************************************* //


std::stack<Particle*> makeStack(Source s , int size) {
  // takes in a pointer to a Source object and the number of histories
  // outputs a particle stack
  
  std::stack<Particle*> pstack; // fite me irl
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
    
//    Particle *p = new Particle(s.location , mu);
    Particle *p = new Particle(0.0 , 1.0);
    // push a pointer to p to the stack
    pstack.push(p);
  }
  return(pstack);
};


// ********************************************************************************************************* //
//  Transport function: Sets up transport problem according to Input object, modifies Tally objects          //
//	                   							         	                                                 //
// ********************************************************************************************************* //

void transport(Leakage &leak_tally , Flux &flux_tally, Slab s , std::stack<Particle*> pstack) 
{
//  this function takes in pointers to a slab and tally objects, and a pointer to a stack of particles, 
//  modifies tally objects tally objects 
//  for each particle in the stack
//  run transport
  double b = s.width / flux_tally.getBins();
  
  int    smallBin;
  int    largeBin;
  double small_x;
  double large_x;

  while( ! pstack.empty() ){  

      while ( pstack.top()->alive ) {
        // do transport
        double optical_dist = -log( Urand() ); // optical distance to collison
        double x1           = pstack.top()->x;

        pstack.top()->x += optical_dist * pstack.top()->mu / s.total_xs;
        double delta_x = pstack.top()->x - x1;
 
        // check for leakage
        if ( pstack.top()->x < 0.0 || pstack.top()->x > s.width ) {
          // leaked out
          pstack.top()->alive = false;
          if ( pstack.top()->x < 0.0 ) {
            // update distance traveled and position to start of the slab
            pstack.top()-> x = 0;
            optical_dist = s.total_xs * x1 / std::abs(pstack.top()->mu); 
          // tally if leaked out of right side
          }
          if ( pstack.top()->x > s.width ) {
            //update leakage tally
            leak_tally.incrementCount();
            // update distance traveled and position to end of the slab 
            pstack.top()->x = s.width;
            optical_dist = s.total_xs * ( pstack.top()->x - x1 ) / std::abs( pstack.top()->mu );
          }
        }
        

        // Flux tally update
        std::vector<double> bounds =  flux_tally.getBounds();
        if( pstack.top()->x >= bounds[0]  and pstack.top()->x <= bounds[1]) {
          // update path flux tally in each spatial bin
          // indexing from 0
          int    oldBin  = floor(x1 / b);
          int    newBin  = floor(pstack.top()->x / b);          
        
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
           
            if (forward == true) {
              smallBin = oldBin;
              largeBin = newBin;
              small_x  = x1;
              large_x  = pstack.top()->x;
            }
            else if (forward == false) {
              smallBin = newBin;
              largeBin = oldBin;
              small_x  = pstack.top()->x;
              large_x  = x1;
            }
          }
          
          // increment the path lengths of the starting and ending mesh bin
          flux_tally.path_lengths[smallBin] += (optical_dist / s.total_xs) * ( ( smallBin +1 ) * b - small_x ) / delta_x;
          flux_tally.path_lengths[largeBin] += (optical_dist / s.total_xs) * (-( largeBin    ) * b + large_x ) / delta_x;
          
          //  increment the path_lengths in all the mesh bins the particle passed through completely
          for ( int i = smallBin + 1 ; i <= largeBin - 1; ++i ) {
            flux_tally.path_lengths[i] += b / std::abs(pstack.top()->mu);
          }
      } // end of flux tally update

          
      // determine next transport step
      if (pstack.top()->alive == true) {
        // have a collision
        if ( Urand() < s.scatter_xs / s.total_xs ) {
          // particle scatters
          pstack.top()->mu = 2.0 * Urand() - 1.0;
        }
        else {
          // absorbed :(
          pstack.top()->alive = false;
        }
      }
    } // end of current particle

    // destruct particle object and remove ptr from stack
    Particle * tmp = pstack.top();
    pstack.pop(); // get rid of the pointer to the particle 
    delete(&tmp); // destruct the memory holding the particle
  } // end of particle stack

};

// ********************************************************************************************************* //
//  Main Function: creates input and tally objects, calls and times transport(), calls output functions      //
//												             //
// ********************************************************************************************************* //

int main()  {
   
  // set up tallies
  int meshBins = 1000;

  // set up transport problem
  Input I; 
  I.nHistories       = 1e7;
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
  s.mu_bins  = { 1.0 };
  s.mu_probs = { 1.0};

  // initialize a particle stack
  std::stack<Particle*> pstack; 
  pstack = makeStack(s , I.nHistories);

  // create and populate a dictionary of Tally objects
  // with the key being the quantity being tallied
  Leakage leak_tally;
  Flux    flux_tally(meshBins , 0  ,I.slab.width);
  
  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  transport(leak_tally, flux_tally, I.slab , pstack);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  
  // finalize tallies
  leak_tally.setVal(I.nHistories);
  flux_tally.setVals(I.nHistories );

  // output results to terminal
  std::cout<< "leakage prob: " << leak_tally.getVal() << std::endl;
  std::cout<< "This run took: " << duration/1000 << " ms" << std::endl;
  
  // output flux to a file
  std::ofstream out;
  out.open("tallies.out");
  std::vector<double> flux = flux_tally.getVals();
  for(int i = 0; i < flux.size(); ++i) {
    out << std::fixed << std::setprecision(7) <<  i* (I.slab.width / flux_tally.getBins() ) << "   ";
    out << std::fixed << std::setprecision(7) << flux[i] << std::endl;
  }

  return(0);
}
