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
    std::vector <double> cellStarts;
    std::vector <double> scatter_xs;
    std::vector <double> total_xs;
    std::vector <double> fission_xs;

    void   setNCells(double ncells);
    double updateCell(Particle *p , double optical_dist);
    double checkLeaks(Particle *p , double optical_dist);
};

void Slab::setNCells(double ncells) {
    nCells = ncells;
};

double Slab::checkLeaks(Particle *p , double delta_x) {
  double delta_x_real = delta_x;
  double xnew = p->x + delta_x;

  if ( xnew  < 0.0 or xnew > width ) {
    // set p.cell to a negative number, to indicate a leak
    p -> cell = -1;

    // leaked out
    p->alive = false;
    if ( xnew < 0.0 ) {
      // update delta_x_real, and position to start of the slab
      delta_x_real = -p->x;
      p->x = 0.0;
    }
    if ( xnew  > width ) {
      // update delta_x_real, and position to end of the slab 
      delta_x_real = width - p->x;
      p->x = width + 0.1;
    }

  }
  return(delta_x_real);
};

double Slab::updateCell(Particle *p , double optical_dist) {
  double delta_x_real = 0;
  double epsilon = 0.00000001;
  double upBound , loBound;

  if(p->cell == cellStarts.size() -1){
    upBound = width;
  }
  else {
      upBound = cellStarts[p->cell+1];
  }
  
  if(p->cell == 0) {
    loBound = 0.0;
  }
  else {
    loBound = cellStarts[p->cell];
  }
  
  double delta_x = optical_dist / total_xs[p->cell] / p->mu;
  
  // check for initial leakage
  delta_x = checkLeaks(p , delta_x);
  if (p-> cell < 0) {
    std::cout << "Leaked! Particle x:" << p->x << " cell:" << p->cell << " delta x: " << delta_x <<std::endl <<std::endl;
    return(delta_x);
  };

  if (p-> x + delta_x < upBound and p->x + delta_x > loBound) {
    p->x += delta_x;
    std::cout << "We didn't change cells!";
    std::cout << " Particle x: " << p->x << " cell: " << p->cell  
              << " delta_x:  "<< delta_x_real << " optical_dist: " << optical_dist  << std::endl;
    return(delta_x);
  };

  // move particle to new cell and update delta_x
  if (p-> x + delta_x > upBound) {
    // if particle went forward a cell
    p->cell+=1; // update cell
    delta_x_real += upBound - p->x; // increment delta_x_real by the distance to the new cell
    p->x = upBound + epsilon; // update particle location
    // recursively call updateCell with optical distance decremented by optical distance traversed 
    // in the current cell
    std::cout << "Moving to the next cell! Particle x: " << p->x << " cell: " << p->cell  
              << " delta_x:  "<< delta_x_real << " optical_dist: " << optical_dist 
              << " OD from move to next cell: " << ((delta_x_real / p->mu) / total_xs[p->cell - 1])
              << std::endl;
    delta_x_real += updateCell(p , optical_dist - ((delta_x_real / p->mu) / total_xs[p->cell - 1]) ); 
  }

  else if (p->x + delta_x < loBound) {
    // if particle went back a celll
    p->cell -= 1;
    delta_x_real -= p->x - loBound; // increment delta_x_real by the distance to the new cell
    p->x = loBound - epsilon; // update particle location
    // recursively call updateCell with optical distance decremented by optical distance traversed 
    // in the current cell
    delta_x_real -= updateCell(p , optical_dist - ((delta_x_real / p->mu) / total_xs[p->cell+1]) );
  }

  return(delta_x_real);

};

// ********************************************************************************************************* //
//  Leakage class, for tallies of leakage probability                                                        //
// ********************************************************************************************************* //

class Leakage {
  public:
    void       setVal(double nHistories);
    void       update(Particle *p , double end);
    double     getVal(void);
               Leakage();

  private:
    double   count;
    double   value;
    double   uncertainty;
};

void Leakage::setVal(double nHistories ) { 
  // Finalizes the leakage tally and calculates uncertainty
  value =  count / nHistories;
  double squared = pow(count , 2);
  uncertainty = sqrt( ( squared / nHistories - pow(value,2) ) / nHistories );
};

void Leakage::update(Particle *p , double end) {
  if(p->x > end and p-> cell < 0) {
    count++;
  }
};


Leakage::Leakage() {
  count = 0;
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
                         Flux(int meshBins , double x0 , double x1);
    std::vector<double>  getBounds(void);
    int                  getBins(void);
    double               getWidth(void);
    std::vector<double>  getVals(void);
    void                 setVals(double nhistories);
    void                 update(Particle *p , double delta_x);
  private:
    double               width;
    int                  nbins;
    std::vector<double>  bounds;
    std::vector <double> path_lengths;
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


void Flux::update(Particle *p , double delta_x) {

     if(p->x >= bounds[0] and p->x <= bounds[1]) {
       // if particle is inside the bounds of the mesh tally
       // update path flux tally in each spatial bin indexing from 0
       double b       = width / nbins;
       int    oldBin  = floor(p->x - delta_x / b);
       int    newBin  = floor(p->x / b);          
       bool   forward = newBin > oldBin; 

       // if the particle stayed in the same bin,
       // add the whole path length to the sum in that bin
       if ( oldBin == newBin ) {
         path_lengths[oldBin] += delta_x / std::abs(p->mu);
       }
   
       else {
         //  otherwise iterate through the bins the particle passed through, 
         //  incrementing them by the path length traversed in that bin
         int    smallBin , largeBin;
         double small_x  , large_x;

         //sort the bins and positions the particles moved between
         if (forward == true) {
           smallBin = oldBin;
           largeBin = newBin;
           small_x  = p->x + delta_x;
           large_x  = p->x;
         }
         else if (forward == false) {
           smallBin = newBin;
           largeBin = oldBin;
           small_x  = p->x;
           large_x  = p->x + delta_x;
         }
       
         // increment the path lengths of the starting and ending mesh bin
         path_lengths[smallBin] +=  ( ( smallBin +1 ) * b - small_x ) / std::abs(p->mu);
         path_lengths[largeBin] +=  (-( largeBin    ) * b + large_x ) / std::abs(p->mu);
       
         //  increment the path_lengths in all the mesh bins the particle passed through completely
         for ( int i = smallBin + 1 ; i <= largeBin - 1; ++i ) {
           path_lengths[i] += b / std::abs(p->mu);
         }
      }
    }

};

// ********************************************************************************************************* //
//  Input class                                                                                              //
// ********************************************************************************************************* //

class Input {
  public:  
    Slab     slab;
    double   nHistories;
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
    Particle *p = new Particle(s.location , mu);
    // push a pointer to p to the stack
    pstack.push(p);
  }

  //std::cout << "Made a particle stack with " << pstack.size() << " particles" << std::endl;
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

  while( ! pstack.empty() ){  
     
      // determine starting cell
      for(int i =0; i < s.cellStarts.size(); ++i) {
        if(pstack.top()->x >= s.cellStarts[i] and pstack.top()->x < s.cellStarts[i+1] ) {
          pstack.top()->cell = i;
        }
      }
      std::cout<<"New Particle!! " <<std::endl;

      while ( pstack.top()->alive ) {
        // do transport

        double optical_dist = -log( Urand() ); // optical distance to collison
     
        std::cout << "Lets take a step" << std::endl;
        std::cout << " x: " << pstack.top()->x <<  " mu: " << pstack.top()->mu << " OD: "<< optical_dist << std::endl; 

        // check if crossed into new cell, and find the total delta_x for the transport
        double delta_x_real = s.updateCell(pstack.top() , optical_dist); 
        
        // tally update
        flux_tally.update(pstack.top() , delta_x_real);
        leak_tally.update(pstack.top() , s.width);
          
        // determine next transport step
        if (pstack.top()->alive == true) {
          // have a collision
          if ( Urand() < s.scatter_xs[pstack.top()->cell] / s.total_xs[pstack.top()->cell] ) {
            // particle scatters
            pstack.top()->mu = 2.0 * Urand() - 1.0;
            std::cout<<"Scatter!!" <<std::endl;
          }
          else {
            // absorbed :(
            std::cout<<"Capture!!" <<std::endl <<std::endl;
            pstack.top()->alive = false;
          }
        }
    } // end of current particle. "He's dead Jim."

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
  I.nHistories       = 1e4;
  I.slab.setNCells(4);
  I.slab.width       = 4.0;
  I.slab.cellStarts  = {0.0 , 1.0  , 2.0 , 3.0};
  I.slab.total_xs    = {1.0 , 1.0  , 1.0 , 1.0};
  I.slab.scatter_xs  = {0.8 , 0.8  , 0.0 , 0.0};
  I.slab.fission_xs  = {0.0 , 0.0  , 0.0 , 0.0};

  // define a source
  Source s;
  s.location =   0.0;
  s.mu_bins  = { 1.0};
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
