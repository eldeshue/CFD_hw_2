// container class for 2D container
// several algorithm for Numerical analysis to solve PDE.
#ifndef _2D_PROFILE_
#define _2D_PROFILE_

#include<vector>
#include<map>
#include<cmath>
#include"Dense_Matrix.h"
#include"Linear_Solver.h"

template< typename T, int time_grid, int space_grid >
class Profile_2D {
  private:

  T start_time;
  T time_step;

  T space_size;
  T space_step;

  T alpha; // coefficient of the term, assumed linear
  T diffusion_factor;

  public:
  std::vector<std::vector<std::pair<double,double>>> data;

  Profile_2D(T start_t, T dt, T h, T dh, T coefficient ) : start_time(start_t), time_step(dt),  
                                                            space_size(h), space_step(dh), alpha(coefficient)
  {
    // set data structure
    data = std::vector<std::vector<std::pair<double,double>>> ( time_grid, std::vector<std::pair<double,double>>( space_grid ) );

    diffusion_factor = alpha * time_step / (space_step * space_step);

  };

  ~Profile_2D(){

  };

  // initializing first time level
	// initial condition, a value
	// t = 0 -> u = init
  void initial_condition(double init){ 
    for( int j = 1; j < space_grid - 1; j++){
		  data[0][j]=std::make_pair(init, j*space_step);
	  }
  };

  // initializing first and last space level
	// boundary condition, two values
	// y_0 -> u = start, y_max -> u = end
  void boundary_condition(double start, double end){
      for( int n = 0; n < time_grid; n++){
		  data[n][0] = std::make_pair(start, 0);
		  data[n][space_grid - 1] = std::make_pair(end,(space_grid - 1)*space_step);
	  }
  };
	// Analytic solution of a simple parabolic equation
  void Para_Analytic(double init_val){
    // calculate velocity profile
	  for( int n = 0; n < time_grid; n++){
		  for( int j = 0; j < space_grid; j++){
        double t = time_step * n;
        double s = space_step * j;
        double eta = s / (2 * sqrt( alpha * t));
        double eta1 = space_size / (2 * sqrt( alpha * t));
        double temp = 0;
        for(int i = 0; i < 100; i++){
          temp += init_val * (std::erfc(2*i*eta1 + eta) - std::erfc(2*(i + 1)*eta1 - eta));
        }
			  data[n][j] = std::make_pair( temp, j * space_step );
		  }
	  }
  }

  // FDE algorithm.
  // Parabolic
  // Forward Time, Central Space Differencing algorithm
  void Para_FTCS_explicit(){
    for( int n = 1; n < time_grid; n++){
		  for( int j = 1; j < space_grid - 1; j++){
			  data[n][j] = std::make_pair(  
						  data[n-1][j].first + diffusion_factor * (data[n-1][j+1].first - 2 * data[n-1][j].first + data[n-1][j-1].first),
						  j * space_step );
		  }
	  }
  }
  // DuFort-Frankel algorithm
  void Para_DF_explicit(){
    for(int j = 1; j < space_grid - 1; j++){
		  data[1][j] = std::make_pair(  
						data[1][j].first + diffusion_factor * (data[1][j+1].first - 2 * data[1][j].first + data[1][j-1].first ) ,
						j * space_step );
	  }

	// n = 2 ~ end
	  for( int n = 2; n < time_grid; n++){
		  for( int j = 1; j < space_grid - 1; j++){
			  data[n][j] = std::make_pair(  
				  ( (1 - 2 * diffusion_factor) * data[n-2][j].first 
				  + 2 * diffusion_factor * (data[n-1][j+1].first + data[n-1][j-1].first)) / (1 + 2 * diffusion_factor) ,
				  j * space_step );
		  }
	  }
  }
  // Laasonen algorithm
	template<int matrix_size>// matrix_size = space_grid - 2;
  void Para_LSN_implicit(){
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
	  for(int n = 1; n < time_grid; n++){
		  // Data from previous Time line, B
		  Dense_Matrix<double, matrix_size, 1> B;
		  // initialize B
		  B.data_table[0][0] = data[n - 1][1].first + diffusion_factor * data[n][0].first;// boundary condition
		  B.data_table[matrix_size - 1][0] = data[n - 1][space_grid - 2].first + diffusion_factor * data[n][space_grid - 1].first;// boundary condition
		  for(int j = 1; j < matrix_size - 1; j++){
			  B.data_table[j][0] = data[n - 1][j + 1].first;
		  }
		  // apply Thomas algorithm and get velocity profiles
		  Dense_Matrix<double, matrix_size, 1> X = Thomas_Algorithm(A, B);
		  for(int j = 1; j < space_grid - 1; j++){
			  data[n][j] = std::make_pair(X.data_table[j - 1][0], j * space_step ); // put each result(velocity) to the data
		  }
	  }

  }
  // Crank - Nicolson algorithm
	template<int matrix_size>// matrix_size = space_grid - 2;
  void Para_CN_implicit(){
		Dense_Matrix<double, matrix_size, matrix_size> A;
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
		for(int n = 1; n < time_grid; n++){
			// Data from previous Time line, B
			Dense_Matrix<double, matrix_size, 1> B;
			// initialize B
			B.data_table[0][0] = // boundary condition
			diffusion_factor * (data[n - 1][2].first - 2 * data[n - 1][1].first + data[n - 1][0].first 
												+ data[n][0].first) + 2 * data[n - 1][1].first;
			B.data_table[matrix_size - 1][0] = // boundary condition
			diffusion_factor * (data[n - 1][space_grid - 1].first - 2 * data[n - 1][space_grid - 2].first + data[n - 1][space_grid - 3].first 
												+ data[n][space_grid - 1].first) + 2 * data[n - 1][space_grid - 2].first;
			for(int j = 1; j < matrix_size - 1; j++){
				B.data_table[j][0] = diffusion_factor * (data[n - 1][j + 2].first - 2 * data[n - 1][j + 1].first + data[n - 1][j].first)
				 + 2 * data[n - 1][j + 1].first;
			}
			// apply Thomas algorithm and get velocity profiles
			Dense_Matrix<double, matrix_size, 1> X = Thomas_Algorithm(A, B);
			for(int j = 1; j < space_grid - 1; j++){
				data[n][j] = std::make_pair(X.data_table[j - 1][0], j * space_step ); // put each result(velocity) to the data
			}
		}
  }
	// after FDE
	void Para_Relative_Error(){
		T init_val = data[0][0].first;
		for(int n = 0; n < time_grid; n++){
			for(int j = 0; j < space_grid; j++){
				T t = time_step * n;
        T s = space_step * j;
        T eta = s / (2 * sqrt( alpha * t));
        T eta1 = space_size / (2 * sqrt( alpha * t));
				T Analytic_value = 0;
        for(int i = 0; i < 100; i++){
          Analytic_value += init_val * (std::erfc(2*i*eta1 + eta) - std::erfc(2*(i + 1)*eta1 - eta));
        }
				data[n][j].first = 100 * (Analytic_value - data[n][j].first) / Analytic_value;
			}
		}
	}

};

#endif