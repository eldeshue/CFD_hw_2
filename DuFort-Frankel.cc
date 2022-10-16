/* 2022.10.14
 * 전산 유체 역학 과제 2
 * solving couette flow, Dufort-Frankel 
 * compile option: g++ -o 파일명 파일명.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "gnuplot-iostream.h" // library for connecting gnuplot
#include "FDE.h"

int main(){
	std::cout << std::setprecision(16);
  // set Gnuplot object;
	Gnuplot gp;

	// set variables
	double U_0 = 40;
	double dynamic_viscous = 0.000217;
	double h = 0.04;
	int JM = 41;
	int NM = 541; // 1) 541    2) 361
	double delta_y = 0.001; 
	double delta_t = 0.002; // 1) 0.002, 2) 0.003 
	double diffusion_factor = dynamic_viscous * delta_t / (delta_y * delta_y);

	// set container for the graph
	// initialize data with initial condition and boundary condition
	// initial condition
  std::vector<std::vector<std::pair<double,double>>> data(NM,std::vector<std::pair<double,double>>(JM));
	int MAX_Y = JM - 1;
	int MAX_N = NM - 1;
	data[0][0] = std::make_pair(U_0, 0);
	for( int j = 1; j < JM; j++){
		data[0][j]=std::make_pair(0, j*delta_y);
	}
	// boundary condition, n = 0
	for( int n = 1; n < NM; n++){
		data[n][0] = std::make_pair(U_0, 0);
		data[n][MAX_Y] = std::make_pair(0, MAX_Y*delta_y);
	}
	// calculate velocity profile with DuFort-Frankel
	// first step FTCS, n = 1
	for(int j = 1; j < MAX_Y; j++){
		data[1][j] = std::make_pair(  
						data[1][j].first + diffusion_factor * (data[1][j+1].first - 2 * data[1][j].first + data[1][j-1].first ) ,
						j * delta_y );
	}
	// n = 2 ~ end
	for( int n = 2; n < NM; n++){
		for( int j = 1; j < MAX_Y; j++){
			data[n][j] = std::make_pair(  
				( (1 - 2 * diffusion_factor) * data[n-2][j].first 
				+ 2 * diffusion_factor * (data[n-1][j+1].first + data[n-1][j-1].first)) / (1 + 2 * diffusion_factor) ,
				j * delta_y );
		}
	}
	
	// make a plot 
	// y value is   and x value is 
	// setting option for the gnuplot
	gp << "set xlabel 'u (m/sec)' font 'Helvetica, 12'\n" << 
    "set ylabel 'y (m)' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title 'Velofity profiles obtained by DuFort-Frankel' font 'Helvetica, 13'" << std::endl;
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
	// NM = 361
	<< gp.file1d(data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(data[60])  << "with linespoints title 't = 0.18'," 
	<< gp.file1d(data[120]) << "with linespoints title 't = 0.36'," 
	<< gp.file1d(data[180]) << "with linespoints title 't = 0.54'," 
	<< gp.file1d(data[240]) << "with linespoints title 't = 0.72'," 
	<< gp.file1d(data[300]) << "with linespoints title 't = 0.90'," 
	<< gp.file1d(data[360]) << "with linespoints title 't = 1.08'"  
	*/
	<< std::endl;
}
