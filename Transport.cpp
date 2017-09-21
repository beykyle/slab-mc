#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Random.h"
#include <string>
#include <algorithm>

using namespace std;
//std::vector<double>

class Particle {
  public:
	double x, mu, x_old, x_new, cell,cell_orig;
    bool   alive;
};

int main() {
  // setup problem
  int    nHistories    = 1e6;
  vector<double>region_boundaries = { 0, 2,3, 4.0 };
  vector<double>total_xs = { 1,8,1};
  vector<double>scatter_ratio = { 0,0.99,0 };
  double thickness     = region_boundaries[size(region_boundaries)-1];
  double mesh_number   = 1000;
  string filename = "test";
 


  double mesh_dist = thickness / mesh_number;  
  vector<double>phi(mesh_number, 0);

  double leakage_hist    = 0.0;
  double leakage         = 0.0;
  double leakage_squared = 0.0;
  Particle p;
  // loop over histories
  for ( int i = 0 ; i < nHistories ; i++ ) { 
    // generate source particle
    p.x     = 0.0;
    p.mu    = 1.0;
	p.cell = 0;
    p.alive = true;

    leakage_hist = 0.0;
    while ( p.alive ) {
      // do transport
      double xsdist_to_collision = -log( Urand() );
	  double xsdist_track = xsdist_to_collision;
	  p.x_old = p.x;
	  p.x += xsdist_to_collision * p.mu/ total_xs[p.cell];
	  if (p.x > (fmax(0, region_boundaries[p.cell])) && p.x < region_boundaries[fmin(size(region_boundaries) - 1,p.cell + 1)]) {
		  p.x = p.x;
	  }
	  else {
		  p.x = p.x_old;
		  if (p.mu < 0) {
			  while (xsdist_track > 0) {

				  xsdist_track -= (p.x - region_boundaries[p.cell]) * total_xs[p.cell];
				  p.x -= (p.x - region_boundaries[p.cell]);
				  p.cell -= 1;
				  if (p.cell < 0) {
					  break;
				  }
			  }
			  p.cell += 1;
			  if (xsdist_track > 0)
				  p.x = -1;
			  else
				  p.x += -1 * xsdist_track / total_xs[p.cell];
		  }
		  else {
			  while (xsdist_track > 0) {
				  p.cell += 1;
				  xsdist_track -= (region_boundaries[p.cell]-p.x) * total_xs[p.cell-1];
				  p.x += (region_boundaries[p.cell]-p.x);
				  if (p.cell == (size(region_boundaries)) - 1) {
					  break;
				  }
			  }
			  p.cell -= 1;
			  if (xsdist_track > 0)
				  p.x = thickness + 1;
			  else
				  p.x -= -1 * xsdist_track / total_xs[p.cell];
		  }
	  }
	  int bin1 = floor(p.x_old / mesh_dist);
	  int bin2 = floor(p.x / mesh_dist);
	  double abs_pmu = abs(p.mu);
	  double full_width_length = mesh_dist / abs_pmu;

	  if (bin1 == bin2) {
		  phi[bin1] += abs(p.x - p.x_old) / p.mu;
	  }
	  else if (bin1 < bin2) {
		  if (p.x > thickness) {
			  phi[bin1] += (mesh_dist*(bin1 + 1) - p.x_old) / abs_pmu;
			  for (int k = bin1 + 1; k < mesh_number; k++) {
				  phi[k] += full_width_length;
			  }
		  }
		  else {
			  phi[bin1] += (mesh_dist*(bin1 + 1) - p.x_old) / abs_pmu;
			  phi[bin2] += (p.x - mesh_dist*bin2) / abs_pmu;
			  for (int k = bin1 + 1; k < bin2; k++) {
				  phi[k] += full_width_length;
			  }
		  }
	  }
	  else if (bin1 > bin2) {
		  if (p.x < 0) {
			  phi[bin1] += (p.x_old - mesh_dist*bin1) / abs_pmu;
			  for (int k = 0; k < bin1; k++) {
				  phi[k] += full_width_length;
			  }
		  }
		  else {
			  phi[bin1] += (p.x_old - mesh_dist*bin1) / abs_pmu;
			  phi[bin2] += (mesh_dist*(bin2 + 1) - p.x) / abs_pmu;
			  for (int k = bin2 + 1; k < bin1; k++) {
				  phi[k] += full_width_length;
			  }
		  }
	  }
      if ( p.x < 0.0 || p.x > thickness ) {
        // leaked out
        p.alive = false;
        // tally if leaked out of right side
        if ( p.x > thickness ) {
          leakage_hist += 1.0;
		}
      }
      else {
		// have a collision
		if (Urand() < scatter_ratio[p.cell]) {
			// particle scatters
			p.mu = 2.0 * Urand() - 1.0;
		}
		else {
			// absorbed :(
			p.alive = false;
		}
      }

    } // end of particle loop

    leakage         += leakage_hist;
    leakage_squared += leakage_hist * leakage_hist;


  } // we are almost done!

  double mean  = leakage / nHistories;
  double stdev = sqrt( ( leakage_squared / nHistories - mean*mean ) / nHistories );

  cout << mean << "  " << stdev / mean << std::endl;

  ofstream phi_file;
  string phi_name = "Phi_" + filename +".txt";
  phi_file.open(phi_name);
  for (int n = 0; n < mesh_number; n++) {
	  phi_file << n*mesh_dist << " " << phi[n]/mesh_dist << endl;
  }
  cout << "Phi File created." << endl;
  phi_file.close();

  system("pause");
  return 0;
}

