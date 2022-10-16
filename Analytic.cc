/* 2022.10.17
 * 전산 유체 역학 과제 2
 * solving couette flow, Analytic solution
 * compile option: g++ -o 파일명 파일명.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include "gnuplot-iostream.h" // library for connecting gnuplot
#include "FDE.h"

const double pi = 3.14159265358979;

int main(){
	std::cout << std::setprecision(16);
  // set Gnuplot object;
	Gnuplot gp;

	// set variables
	double U_0 = 40;
	double dynamic_viscous = 0.000217;
	double h = 0.04;
	int JM = 41;
	int NM = 541; // 1) 541    2) 541
	double delta_y = 0.001; 
	double delta_t = 0.002; // 1) 0.002, 2) 0.00232 
	double diffusion_factor = dynamic_viscous * delta_t / (delta_y * delta_y);

	// set container for the graph
	// initialize data with initial condition and boundary condition
	// initial condition
  std::vector<std::vector<std::pair<double,double>>> data(NM,std::vector<std::pair<double,double>>(JM));
	int MAX_Y = JM - 1;
	int MAX_N = NM - 1;
	
	// calculate velocity profile
	for( int n = 0; n < NM; n++){
		for( int j = 0; j < MAX_Y; j++){
      double t = delta_t * n;
      double y = delta_y * j;
      double eta = y/ (2 * sqrt( dynamic_viscous * t));
      double eta1 = h / (2 * sqrt( dynamic_viscous * t));
      double temp = 0;
      for(int i = 0; i < 100; i++){
        temp += U_0 * (std::erfc(2*i*eta1 + eta) - std::erfc(2*(i + 1)*eta1 - eta));
      }
			data[n][j] = std::make_pair( temp, j * delta_y );
		}
	}
	// make a plot
	// find row number of data, and print that row. 
	// y value is   and x value is 
	// setting option for the gnuplot
	gp << "set xlabel 'u (m/sec)' font 'Helvetica, 12'\n" << 
    "set ylabel 'y (m)' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title 'Velofity profiles obtained by Analytic' font 'Helvetica, 13'" << std::endl;
	gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" << 
    "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
	gp << "set autoscale \n"; 
  gp << "set grid \n";
  gp << "plot" 
	<< gp.file1d(data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(data[90])  << "with linespoints title 't = 0.18'," 
	<< gp.file1d(data[180]) << "with linespoints title 't = 0.36',"
	<< gp.file1d(data[270]) << "with linespoints title 't = 0.54',"
	<< gp.file1d(data[360]) << "with linespoints title 't = 0.72',"
	<< gp.file1d(data[480]) << "with linespoints title 't = 0.90',"
	<< gp.file1d(data[540]) << "with linespoints title 't = 1.08'" 
	/*
	<< gp.file1d(data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(data[90])  << "with linespoints title 't = 0.21',"
	<< gp.file1d(data[180]) << "with linespoints title 't = 0.42',"
	<< gp.file1d(data[270]) << "with linespoints title 't = 0.63',"
	<< gp.file1d(data[360]) << "with linespoints title 't = 0.84',"
	<< gp.file1d(data[480]) << "with linespoints title 't = 1.04',"
	<< gp.file1d(data[540]) << "with linespoints title 't = 1.25'"
	*/
	<< std::endl;
}
