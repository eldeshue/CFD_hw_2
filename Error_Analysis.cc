// 2022.10.16 CFD hw2 
// Error Analysis for solving Parabolic eq, FDE
// compile option: g++ -o Error_Analysis Error_Analysis.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include "gnuplot-iostream.h" // library for connecting gnuplot
#include "FDE.h"
#include "Dense_Matrix.h"
#include "Linear_Solver.h"

int main(){

  // set variables
	double U_0 = 40;
	double dynamic_viscous = 0.000217;
	double h = 0.04;
	int JM = 41;
	int NM = 541; // 91 step -> 0.18 sec, 541 step -> 1.08 sec.
	double delta_y = 0.001; 
	double delta_t = 0.002;  
	double diffusion_factor = dynamic_viscous * delta_t / (delta_y * delta_y);
  int MAX_Y = JM - 1;
	int MAX_N = NM - 1;

  // calculating Exact data through complimentary Errorfunction, erfc
  // 2D array
  double Analytic_value[NM][JM];
  for( int n = 0; n < NM; n++ ){
    for( int j = 0; j < JM; j++ ){
      double t = delta_t * n;
      double y = delta_y * j;
      double eta = y/ (2 * sqrt( dynamic_viscous * t));
      double eta1 = h / (2 * sqrt( dynamic_viscous * t));
      double temp = 0;
      for(int i = 0; i < 100; i++){
        temp += U_0 * (std::erfc(2*i*eta1 + eta) - std::erfc(2*(i + 1)*eta1 - eta));
      }
      Analytic_value[n][j] = temp;
    }    
  }

  // claculating FDE
  // calculating FDE through FTCS
  std::vector<std::vector<std::pair<double,double>>> FTCS(NM,std::vector<std::pair<double,double>>(JM));
	FTCS[0][0] = std::make_pair(U_0, 0);
	for( int j = 1; j < JM; j++){
		FTCS[0][j]=std::make_pair(0, j*delta_y);
	}
	// boundary condition, n = 0
	for( int n = 1; n < NM; n++){
		FTCS[n][0] = std::make_pair(U_0, 0);
		FTCS[n][MAX_Y] = std::make_pair(0, MAX_Y*delta_y);
	}

  for( int n = 1; n < NM; n++){
		for( int j = 1; j < MAX_Y; j++){
			FTCS[n][j] = std::make_pair(  
						FTCS[n-1][j].first + diffusion_factor * (FTCS[n-1][j+1].first - 2 * FTCS[n-1][j].first + FTCS[n-1][j-1].first),
						j * delta_y );
		}
	}
    
  // calculating FDE through DuFort-Frankel
  std::vector<std::vector<std::pair<double,double>>> DF(NM,std::vector<std::pair<double,double>>(JM));
	DF[0][0] = std::make_pair(U_0, 0);
	for( int j = 1; j < JM; j++){
		DF[0][j]=std::make_pair(0, j*delta_y);
	}
	// boundary condition, n = 0
	for( int n = 1; n < NM; n++){
		DF[n][0] = std::make_pair(U_0, 0);
		DF[n][MAX_Y] = std::make_pair(0, MAX_Y*delta_y);
	} 
  // calculate velocity profile with DuFort-Frankel
	// first step FTCS, n = 1
	for(int j = 1; j < MAX_Y; j++){
		DF[1][j] = std::make_pair(  
						DF[1][j].first + diffusion_factor * (DF[1][j+1].first - 2 * DF[1][j].first + DF[1][j-1].first ) ,
						j * delta_y );
	}
	// n = 2 ~ end
	for( int n = 2; n < NM; n++){
		for( int j = 1; j < MAX_Y; j++){
			DF[n][j] = std::make_pair(  
				( (1 - 2 * diffusion_factor) * DF[n-2][j].first 
				+ 2 * diffusion_factor * (DF[n-1][j+1].first + DF[n-1][j-1].first)) / (1 + 2 * diffusion_factor) ,
				j * delta_y );
		}
	}
  
  // calculating FDE through Laasonen
  std::vector<std::vector<std::pair<double,double>>> LSN(NM,std::vector<std::pair<double,double>>(JM));
	LSN[0][0] = std::make_pair(U_0, 0);
	for( int j = 1; j < JM; j++){
		LSN[0][j]=std::make_pair(0, j*delta_y);
	}
	// boundary condition, n = 0
	for( int n = 1; n < NM; n++){
		LSN[n][0] = std::make_pair(U_0, 0);
		LSN[n][MAX_Y] = std::make_pair(0, MAX_Y*delta_y);
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
		// LSN from previous Time line, B
		Dense_Matrix<double, matrix_size, 1> B;
		// initialize B
		B.data_table[0][0] = LSN[n - 1][1].first + diffusion_factor * LSN[n][0].first;// boundary condition
		B.data_table[matrix_size - 1][0] = LSN[n - 1][MAX_Y - 1].first + diffusion_factor * LSN[n][MAX_Y].first;// boundary condition
		for(int j = 1; j < matrix_size - 1; j++){
			B.data_table[j][0] = LSN[n - 1][j + 1].first;
		}
		// apply Thomas algorithm and get velocity profiles
		Dense_Matrix<double, matrix_size, 1> X = Thomas_Algorithm(A, B);
		for(int j = 1; j < MAX_Y; j++){
			LSN[n][j] = std::make_pair(X.data_table[j - 1][0], j * delta_y ); // put each result(velocity) to the LSN
		}
	}
  
  // calculating FDE through Crank-Nicolson.
  std::vector<std::vector<std::pair<double,double>>> CN(NM,std::vector<std::pair<double,double>>(JM));
	CN[0][0] = std::make_pair(U_0, 0);
	for( int j = 1; j < JM; j++){
		CN[0][j]=std::make_pair(0, j*delta_y);
	}
	// boundary condition, n = 0
	for( int n = 1; n < NM; n++){
		CN[n][0] = std::make_pair(U_0, 0);
		CN[n][MAX_Y] = std::make_pair(0, MAX_Y*delta_y);
	} 
  for(int i = 0; i < matrix_size; i++ ){
		for(int j = 0; j < matrix_size; j++){
			if(j == i){
				A.data_table[i][j] = 2 * (1 + diffusion_factor);
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
		// CN from previous Time line, B
		Dense_Matrix<double, matrix_size, 1> B;
		// initialize B
		B.data_table[0][0] = // boundary condition
		diffusion_factor * (CN[n - 1][2].first - 2 * CN[n - 1][1].first + CN[n - 1][0].first 
											+ CN[n][0].first) + 2 * CN[n - 1][1].first;
		B.data_table[matrix_size - 1][0] = // boundary condition
		diffusion_factor * (CN[n - 1][MAX_Y].first - 2 * CN[n - 1][MAX_Y - 1].first + CN[n - 1][MAX_Y - 2].first 
											+ CN[n][MAX_Y].first) + 2 * CN[n - 1][MAX_Y - 1].first;
		for(int j = 1; j < matrix_size - 1; j++){
			B.data_table[j][0] = diffusion_factor * (CN[n - 1][j + 2].first - 2 * CN[n - 1][j + 1].first + CN[n - 1][j].first)
			 + 2 * CN[n - 1][j + 1].first;
		}
		// apply Thomas algorithm and get velocity profiles
		Dense_Matrix<double, matrix_size, 1> X = Thomas_Algorithm(A, B);
		for(int j = 1; j < MAX_Y; j++){
			CN[n][j] = std::make_pair(X.data_table[j - 1][0], j * delta_y ); // put each result(velocity) to the CN
		}
	}
  // Error Calculation
  for( int n = 0; n < NM; n++){
    for(int j = 0; j < JM; j++){
      
    }
  }
  
  // Error Calculation
  for( int n = 0; n < NM; n++){
    for(int j = 0; j < JM; j++){
      FTCS[n][j].first = 100 * ( Analytic_value[n][j] - FTCS[n][j].first) / Analytic_value[n][j];
      DF[n][j].first = 100 * ( Analytic_value[n][j] - DF[n][j].first) / Analytic_value[n][j];
      LSN[n][j].first = 100 * ( Analytic_value[n][j] - LSN[n][j].first) / Analytic_value[n][j];
      CN[n][j].first = 100 * ( Analytic_value[n][j] - CN[n][j].first) / Analytic_value[n][j];
    }
  }
  
  // print these data through Gnuplot.
  // set Gnuplot object;
	Gnuplot gp;
  gp << "set xlabel 'Error (%)' font 'Helvetica, 12'\n" << 
    "set ylabel 'y (m)' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title ' Comparison of error distributions, t = 0.18 sec ' font 'Helvetica, 13'" << std::endl;
	gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" << 
    "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
	gp << "set autoscale \n"; 
  gp << "set grid \n";
  gp << "plot" 
	// NM = 541, delta_t = 0.002
	<< gp.file1d(FTCS[540])   << "with linespoints title 'FTCS',"
  << gp.file1d(DF[540])  << "with linespoints title 'D/F'," 
	<< gp.file1d(LSN[540]) << "with linespoints title 'LAASONEN'," 
	<< gp.file1d(CN[540]) << "with linespoints title 'C/N'" 
  << std::endl;

}