/* 2022.10.17
 * 전산 유체 역학 과제 2
 * solving couette flow, Error Analysis by step size
 * compile option1: simple error
 g++ -D_ERROR_ -o Error_Analysis_Laasonen Error_Analysis_Laasonen.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 * compile option2: relative error
 g++ -D_RELATIVE_ERROR_ -o Error_Analysis_Laasonen Error_Analysis_Laasonen.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */
 

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include "2D_Profile.h"
#include "gnuplot-iostream.h" // library for connecting gnuplot
#include "FDE.h"

int main(){
  // set variables
	double U_0 = 40;
	double dynamic_viscous = 0.000217;
	double h = 0.04;
	const int JM = 41;
	const int NM1 = 1001;
  const int NM2 = 201;
  const int NM3 = 11;
  const int NM4 = 6; 
	double delta_y = 0.001; 
	double delta_t1 = 0.005;
  double delta_t2 = 0.01;
  double delta_t3 = 0.1;
  double delta_t4 = 0.2;  

  Profile_2D<double, NM1, JM> test1(0, delta_t1, h, delta_y, dynamic_viscous);
  test1.initial_condition(0);
  test1.boundary_condition(U_0, 0);
  test1.Para_LSN_implicit<39>();
  
  Profile_2D<double, NM2, JM> test2(0, delta_t2, h, delta_y, dynamic_viscous);
  test2.initial_condition(0);
  test2.boundary_condition(U_0, 0);
  test2.Para_LSN_implicit<39>();

  Profile_2D<double, NM3, JM> test3(0, delta_t3, h, delta_y, dynamic_viscous);
  test3.initial_condition(0);
  test3.boundary_condition(U_0, 0);
  test3.Para_LSN_implicit<39>();

  Profile_2D<double, NM4, JM> test4(0, delta_t4, h, delta_y, dynamic_viscous);
  test4.initial_condition(0);
  test4.boundary_condition(U_0, 0);
  test4.Para_LSN_implicit<39>();

  // Error Analysis
  // simple error
  #ifdef _ERROR_
  test1.Para_Error();
  test2.Para_Error();
  test3.Para_Error();
  test4.Para_Error();
  #endif
  #ifdef _RELATIVE_ERROR_
  // relative error, normalize error
  test1.Para_Relative_Error();
  test2.Para_Relative_Error();
  test3.Para_Relative_Error();
  test4.Para_Relative_Error();
  #endif
  // print these data through Gnuplot.
  // set Gnuplot object;
	Gnuplot gp;
  gp << "set xlabel 'Error' font 'Helvetica, 12'\n" << 
    "set ylabel 'y (m)' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title ' Effect of time steps, t = 1 sec, Laasonen ' font 'Helvetica, 13'" << std::endl;
	gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" << 
    "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
	gp << "set autoscale \n"; 
  gp << "set grid \n";
  gp << "plot" 
	// NM = 541, delta_t = 0.002
	<< gp.file1d(test1.data[NM1 - 1]) << "with linespoints title 'dt = 0.005',"
  << gp.file1d(test2.data[NM2 - 1]) << "with linespoints title 'dt = 0.01'," 
	<< gp.file1d(test3.data[NM3 - 1]) << "with linespoints title 'dt = 0.1'," 
	<< gp.file1d(test4.data[NM4 - 1]) << "with linespoints title 'dt = 0.2'" 
  << std::endl;

}