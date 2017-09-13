#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <numeric>
//comment
//TODO finish Input.cpp
//TODO make slab.cpp
#include "Random.h"

class Slab {
  public:
    double               width;
    int                  nBins;
    std::vector <double> scatter_xs;
    std::vector <double> total_xs;
    std::vector <double> fission_xs;
};

class Input {
  public:  
    Slab    slab;
    double  nHistories;
};

class Tally {
  public:
    bool                 vector;
    bool                 total;
    double               value , uncertainty;
    std::vector <double> values, uncertainties;
};


class Particle {
  public:
    double x, mu;
    bool   alive;
};

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
  
  double b = slab.width / slab.nBins;

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

std::map<std::string , Tally> transport(Input I) {
  
  // initialaize variables from input object
  int    nHistories    = I.nHistories;
  double total_xs      = I.slab.total_xs[0];
  double scatter_ratio = I.slab.scatter_xs[0] / total_xs;
  double thickness     = I.slab.width;
  
  double nBins           = I.slab.nBins;
  double b               = thickness / nBins;

  
  // initialize counters for leakage tally
  double leakage_hist    = 0.0;
  double leakage         = 0.0;
  double leakage_squared = 0.0;

  // initialize counters for flux tally
  std::vector <double> path_lengths(nBins , 0.0);
  std::vector <double> path_lengths_squared(nBins , 0.0);

  // loop over histories
  for ( int i = 0 ; i < I.nHistories ; i++ ) { 
    // generate source particle
    Particle p;
    p.x     = 0.0;
    p.mu    = 1.0;
    p.alive = true;

    leakage_hist = 0.0;

    while ( p.alive ) {
      // do transport

      double dist_to_collision = -log( Urand() ) / total_xs;
      double x1                = p.x;

      p.x += dist_to_collision * p.mu;
      double delta_x           = p.x - x1;

      // check for leakage
      if ( p.x < 0.0 || p.x > thickness ) {
        // leaked out
        p.alive = false;
        if ( p.x > thickness ) {
          //update leakage tally
          leakage_hist += 1.0;
        }
      }

      
      // determine next transport step
      if (p.alive == true) {
        // have a collision
        if ( Urand() < scatter_ratio ) {
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
    leakage         += leakage_hist;
    leakage_squared += leakage_hist * leakage_hist;


  } // we are almost done!

  
  // calculate tallies and create tally objects
  // leakage
  Tally leak_tally;
  leak_tally.value =  leakage / nHistories;
  leak_tally.uncertainty = sqrt( ( leakage_squared / nHistories - pow(leak_tally.value,2) ) / nHistories );
  // flux
  Tally flux_tally;
  double  sum_path_lengths = 0;
  std::vector <double> diff_path_lengths_sqrd(nBins , 0);

  for(int i = 0; i < path_lengths.size(); ++i) {
    path_lengths_squared[i] = path_lengths[i] * path_lengths[i];
    sum_path_lengths += path_lengths[i];
    flux_tally.values.push_back(path_lengths[i] / b / nHistories);
    flux_tally.uncertainties.push_back(0); 
  }

  flux_tally.value = sum_path_lengths;
  flux_tally.uncertainty = 0; 
  
  // set flags for output formating 
  leak_tally.vector = false;
  leak_tally.total  = false;
  flux_tally.vector = true;
  flux_tally.total  = true;
  
  //  create and populate a dictionary of tallies
  std::map<std::string , Tally> tallies;
  tallies["Leakage Probability"] = leak_tally;
  tallies["Flux"] = flux_tally;
  return(tallies);
};

int main()  {
  
  // set up transport problem
  Input I; 
  I.nHistories       = 1e7;
  I.slab.total_xs    = {1.0};
  I.slab.scatter_xs  = {0.9};
  I.slab.fission_xs  = {0.0};
  I.slab.width       = 4.0;
  I.slab.nBins       = 1000;
   
  // create a dictionary holding the tallies, with the key being the quantity being tallied
  std::map<std::string , Tally>  tallies;
  
  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  tallies = transport(I);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  
  // output results 
  terminal_out(tallies , duration);
  file_out(I.slab , tallies , duration); 
  
  return(0);
}
