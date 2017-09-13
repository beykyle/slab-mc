#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <numeric>

#include "Random.h"

//TODO make input class  in Input.cpp and Input.h, with read from file function

class Slab {
  public:
    double               width;
    int                  nCells;
    std::vector <double> scatter_xs;
    std::vector <double> total_xs;
    std::vector <double> fission_xs;
};

class Input {
  public:  
    Slab    slab;
    int     tallyBins;
    double  nHistories;
};

class Leakage {
  public:
    double leakage_squared;
    double leakahe_hist;
    //function to implement leakage counting
};

class Flux {
  public:
    //functions to implement flux_tally
};

class Tally {
  public:
    bool                 vector;
    bool                 total;
    double               value , uncertainty;
    std::vector <double> values, uncertainties;
    int nbins;
};


class Particle {
  public:
    double x, mu;
    int    cell;
    bool   alive;
};

// ********************************************************************************************************* //
//  Output functions: take in a dictionary of Tally objects, and the duration of time required to run        //
//      transport(). Outputs results to terminal and an output file "tallies.out"                            //
// ********************************************************************************************************* //

void terminal_out(std::map<std::string , Tally> tallies , double duration) {
// a function to output the results of the simulation to the terminal
// requires: a dictionary of Tally objects

  //  iterate over tally dictionary, output tally results to terminal
  std::cout << std::endl  <<"============================================"  << std::endl;
   
  for (const auto &keyVal : tallies ) {
    if (keyVal.second.vector == false or keyVal.second.total == true) {
      std::string total = (keyVal.second.total == true) ? ("Total "):("");
      std::cout << std::setw (25 - keyVal.first.length()) << total << keyVal.first << ": " ;
      std::cout << std::fixed << std::setprecision(6) << keyVal.second.value; 
      std::cout << std::setw (4) << "   stdev: " << std::fixed << std::setprecision(6);
      std::cout << keyVal.second.uncertainty << std::endl;
    }
  } 
  
  std::cout << std::endl  <<"============================================"  << std::endl;
  
  //  output timing results to terminal
  std::cout << "This run took "  << duration / 1000 << " miliseconds"<< std::endl;
  
  std::cout << "============================================"  << std::endl;
  
};

void file_out(Slab slab , std::map<std::string , Tally> tallies , double duration) {
// a function to output the vector tallies  of the simulation to a file
// requires: a dictionary of Tally objects
  
  std::ofstream out;
  out.open("tallies.out");
  out << std::endl  <<"============================================"  << std::endl;
  out << "This run took "  << duration / 1000 << " miliseconds"<< std::endl; 
  out << std::endl  <<"============================================"  << std::endl;

  //  iterate over tally dictionary, output tally results to terminal
  for (const auto &keyVal : tallies ) {
    if (keyVal.second.vector == false or keyVal.second.total == true) {
      std::string total = (keyVal.second.total == true) ? ("Total "):("");
      out  << total << keyVal.first << ": " << keyVal.second.value;
      out  << "    stdev: " << keyVal.second.uncertainty << std::endl;
      // output vector results
      if (keyVal.second.vector == true) {
        out << "x       Value     Uncertainty" << std::endl;
        double b = slab.width / keyVal.second.nbins;
        for (int i = 0; i < keyVal.second.values.size(); ++i) {
          out << std::fixed << std::setprecision(7) << b * i << "   ";
          out << std::fixed << std::setprecision(7) << keyVal.second.values[i];
	  out << "      " << keyVal.second.uncertainties[i] << std::endl;
        }
      }; // loop over tally values
      out << std::endl  <<"============================================"  << std::endl;
    }
  } // loop over tallies
};

// ********************************************************************************************************* //
//  Transport function: Sets up transport problem according to Input object, modifies Tally objects          //
//	                   							         	             //
// ********************************************************************************************************* //

void transport(Tally &leak_tally , Tally &flux_tally, Input I) {
// this function takes in an input object, runs transport, and outputs a dictionary of tally objects 
  
  //// initialize counters for tally objects ////
  // leakage
  double leakage_squared  = 0.0;
  leak_tally.value        = 0.0;
  leak_tally.uncertainty  = 0.0;
  // flux
  double               sum_path_lengths = 0;
  std::vector <double> path_lengths(flux_tally.nbins , 0.0);
  std::vector <double> path_lengths_squared(flux_tally.nbins , 0.0);

  // set flags on tally objects for output formating 
  leak_tally.vector = false;
  leak_tally.total  = false;
  flux_tally.vector = true;
  flux_tally.total  = true;
  
  //// initialaize variables from input object ////
  int    nHistories      = I.nHistories;
  
  double thickness       = I.slab.width;
  double cellWidth       = thickness / I.slab.nCells;       // thickness of bins for xsec sampling
  double tWidth          = thickness / I.tallyBins;  // tickness of mesh bins for tallies

  std::vector <double> scatter_ratio(I.slab.nCells , 0.0);
  for ( int i = 0; i < I.slab.nCells; ++i) {
    scatter_ratio[i] = I.slab.scatter_xs[i] / I.slab.total_xs[i];
  }
 
  //// loop over histories ////
  for ( int i = 0 ; i < I.nHistories ; i++ ) { 
    // generate source particle
    Particle p;
    p.x     = 0.0;
    p.mu    = 1 - Urand(); //sample an isotropic incident particle in 2pi
    p.cell  = 0;
    p.alive = true;
    
    while ( p.alive ) {
      // do transport for a single particle
      
      // sample new distance to collision
      double dist_to_collision = -log( Urand() ) / I.slab.total_xs[p.cell];
      // update particle location and cell
      p.x += dist_to_collision * p.mu;
      p.cell = floor(p.x / cellWidth);

      // check for leakage
      if ( p.x < 0.0 || p.x > thickness ) {
        // leaked out
        p.alive = false;
        if ( p.x > thickness ) {
          //update leakage tally
          leak_tally.value += 1.0;
        }
      }
      
      // determine next transport step
      if (p.alive == true) {
        // have a collision
        if ( Urand() < scatter_ratio[p.cell] ) {
          // particle scatters
          p.mu = 2.0 * Urand() - 1.0;
        }
        else {
          // absorbed :(
          p.alive = false;
        }
      }
    } // end of particle loop
    
    //  update tallies
    leakage_squared += pow(leak_tally.value,2);

  } // we are almost done!
  
  //// calculate tallies ////
  
  // leakage
  leak_tally.value =  leak_tally.value / nHistories;
  leak_tally.uncertainty = sqrt( ( leakage_squared / nHistories - pow(leak_tally.value,2) ) / nHistories );
  
  // flux
  for(int i = 0; i < path_lengths.size(); ++i) {
    path_lengths_squared[i] = path_lengths[i] * path_lengths[i];
    sum_path_lengths += path_lengths[i];
    flux_tally.values.push_back(path_lengths[i] / tWidth / nHistories);
    flux_tally.uncertainties.push_back(0); 
  }

  flux_tally.value = sum_path_lengths;
  flux_tally.uncertainty = 0; 
};

// ********************************************************************************************************* //
//  Main Function: creates input and tally objects, calls and times transport(), calls output functions      //
//												             //
// ********************************************************************************************************* //

int main()  {
  
  // set up transport problem
  Input I; 
  I.nHistories       = 1e7;
  I.slab.total_xs    = {1.0 , 1.0 , 1.0, 1.0 , 1.0};
  I.slab.scatter_xs  = {0.0 , 0.0 , 0.8 , 0.8 , 0.5};
  I.slab.fission_xs  = {0.0};
  I.slab.width       = 10.0;
  I.slab.nCells      = 5;
  
  // Initialize tally objects 
  Tally leak_tally;
  leak_tally.nbins = 1;
  
  Tally flux_tally;
  flux_tally.nbins = 1000;
  
  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  transport(leak_tally, flux_tally, I);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  // create and populate a dictionary of Tally objects
  // with the key being the quantity being tallied
  std::map<std::string , Tally> tallies;
  tallies["Leakage Probability"] = leak_tally;
  tallies["Flux"] = flux_tally;
  
  // output results 
  terminal_out(tallies , duration);
  file_out(I.slab , tallies , duration); 
  
  return(0);
}
