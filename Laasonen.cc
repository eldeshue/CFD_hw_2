/* 2022.10.14
 * 전산 유체 역학 과제 2
 * solving couette flow, Laasonen
 * compile option: g++ -o 파일명 파일명.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "gnuplot-iostream.h" // library for connecting gnuplot
#include "FDE.h"
#include "Dense_Matrix.h"
#include "Linear_Solver.h"

int main(){
	std::cout << std::setprecision(16);
  // set Gnuplot object;
	Gnuplot gp;

	// set variables
	double U_0 = 40;
	double dynamic_viscous = 0.000217;
	double h = 0.04;
	int JM = 41;
	int NM = 541; // 1) 541    2) 109
	double delta_y = 0.001; 
	double delta_t = 0.002; // 1) 0.002, 2) 0.1 
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
	
	const int matrix_size = 39;// Matrix size = JM(41) - 2 (boundary condition) = 39;
	// translate linear algebraic eq to matrix form, AX = B
	// Coefficient Matrix A, tridiagonal
	Dense_Matrix<double, matrix_size, matrix_size> A;
	for(int i = 0; i < matrix_size; i++ ){
		for(int j = 0; j < matrix_size; j++){
			if(j == i){
				A.data_table[i][j] = 1 + 2 * diffusion_factor;
			}
			else if(j == (i + 1) || j == (i - 1)){
				A.data_table[i][j] = 0 - diffusion_factor;
			}
			else{
				A.data_table[i][j] = 0;
			}
		}
	}

	// Calculate Velocity Profile
	for(int n = 1; n < NM; n++){
		// Data from previous Time line, B
		Dense_Matrix<double, matrix_size, 1> B;
		// initialize B
		B.data_table[0][0] = data[n - 1][1].first + diffusion_factor * data[n][0].first;// boundary condition
		B.data_table[matrix_size - 1][0] = data[n - 1][MAX_Y - 1].first + diffusion_factor * data[n][MAX_Y].first;// boundary condition
		for(int j = 1; j < matrix_size - 1; j++){
			B.data_table[j][0] = data[n - 1][j + 1].first;
		}
		// apply Thomas algorithm and get velocity profiles
		Dense_Matrix<double, matrix_size, 1> X = Thomas_Algorithm(A, B);
		for(int j = 1; j < MAX_Y; j++){
			data[n][j] = std::make_pair(X.data_table[j - 1][0], j * delta_y ); // put each result(velocity) to the data
		}
	}
	
	// make a plot
	// find row number of data, and print that row. 
	// y value is   and x value is 
	// setting option for the gnuplot
	gp << "set xlabel 'u (m/sec)' font 'Helvetica, 12'\n" << 
    "set ylabel 'y (m)' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title 'Velofity profiles obtained by Laasonen' font 'Helvetica, 13'" << std::endl;
	gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" << 
    "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
	gp << "set autoscale \n"; 
  gp << "set grid \n";
  gp << "plot" 
	// NM = 541, delta_t = 0.002
	<< gp.file1d(data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(data[90])  << "with linespoints title 't = 0.18'," 
	<< gp.file1d(data[180]) << "with linespoints title 't = 0.36'," 
	<< gp.file1d(data[270]) << "with linespoints title 't = 0.54'," 
	<< gp.file1d(data[360]) << "with linespoints title 't = 0.72'," 
	<< gp.file1d(data[480]) << "with linespoints title 't = 0.90'," 
	<< gp.file1d(data[540]) << "with linespoints title 't = 1.08'" 
	/*
	// NM = 109, delta_t = 0.01
	<< gp.file1d(data[0])   << "with linespoints title 't = 0.00',"
  << gp.file1d(data[18])  << "with linespoints title 't = 0.18',"
	<< gp.file1d(data[36]) << "with linespoints title 't = 0.36',"
	<< gp.file1d(data[54]) << "with linespoints title 't = 0.54',"
	<< gp.file1d(data[72]) << "with linespoints title 't = 0.72',"
	<< gp.file1d(data[90]) << "with linespoints title 't = 0.90',"
	<< gp.file1d(data[108]) << "with linespoints title 't = 1.08'"
	*/
	<< std::endl;
}
