/* 2022.10.17
 * 전산 유체 역학 과제 2
 * solving couette flow, Analytic solution
 * compile option: g++ -o test_profile test_profile.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include "2D_Profile.h"
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
	const int JM = 41;
	const int NM = 541; // 1) 541    2) 541
	double delta_y = 0.001; 
	double delta_t = 0.002; // 1) 0.002, 2) 0.00232 
	
  Profile_2D<double, NM, JM > test(0, delta_t, h, delta_y, dynamic_viscous);
  //test.Para_Analytic(U_0);
  test.initial_condition(0);
  test.boundary_condition(U_0, 0);
  //test.Para_FTCS_explicit();
  //test.Para_DF_explicit();
  test.Para_LSN_implicit();
  
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
	<< gp.file1d(test.data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(test.data[90])  << "with linespoints title 't = 0.18'," 
	<< gp.file1d(test.data[180]) << "with linespoints title 't = 0.36',"
	<< gp.file1d(test.data[270]) << "with linespoints title 't = 0.54',"
	<< gp.file1d(test.data[360]) << "with linespoints title 't = 0.72',"
	<< gp.file1d(test.data[480]) << "with linespoints title 't = 0.90',"
	<< gp.file1d(test.data[540]) << "with linespoints title 't = 1.08'" 
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
