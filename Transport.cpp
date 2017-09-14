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
#include "boost/variant.hpp"

//TODO make input class  in Input.cpp and Input.h, with read from file function

class Particle {
  public:
    double x, mu;
    int    cell;
    bool   alive;
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
    Slab                          slab;
    int                           tallyBins;
    double                        nHistories;
    std::map<std::string , boost::variant <Leakage , Flux > > tallies;
};


// ********************************************************************************************************* //
//  Output functions: take in a dictionary of Tally objects, and the duration of time required to run        //
//      transport(). Outputs results to terminal and an output file "tallies.out"                            //
// ********************************************************************************************************* //

void terminal_out(std::map<std::string , boost::variant <Leakage , Flux > > tallies , double duration) {
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

void file_out(Slab slab , std::map<std::string , boost::variant <Leakage , Flux > > tallies , double duration) {
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

void transport(Tally &leak_tally , Tally &flux_tally, Slab s) {
// this function takes in a slab object, a number of histories, runs transport, and outputs a dictionary of tally objects 
  
};

// ********************************************************************************************************* //
//  Main Function: creates input and tally objects, calls and times transport(), calls output functions      //
//												             //
// ********************************************************************************************************* //

int main()  {
  
  // set up transport problem
  //TODO create input class constructor that parses from input file
  Input I; 
  I.nHistories       = 1e7;
  I.slab.width       = 4.0;
  I.slab.nCells      = 4;
  I.slab.total_xs    = {0.0 , 0.3 , 1.0 , 1.0 } ;
  I.slab.scatter_xs  = {0.0 , 0.25 , 0.7 , 0.1};
  I.slab.fission_xs  = {0.0};

  // create and populate a dictionary of Tally objects
  // with the key being the quantity being tallied
  Leakage leak_tally;
  Flux    flux_tally;
  flux_tally.nbins = 1000; //the mesh for the flux tally has 1000 bins

  // set flags for output formating 
  leak_tally.vector = false;
  leak_tally.total  = false;
  flux_tally.vector = true;
  flux_tally.total  = true;

  I.tallies["Leakage Probability"] = leak_tally;
  I.tallies["Flux"] = flux_tally;
  
  // run and time the transport function
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  transport(I.tallies["Leakage Probability"], I.tallies["Flux"], I.slab);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  
  // output results 
  terminal_out( I.tallies , duration);
  file_out(I.slab , I.tallies , duration); 
  
  return(0);
}
