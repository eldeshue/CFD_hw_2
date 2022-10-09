/*
* 2022.10.29 
* Advanced Matrix Library, made by Dohoon Gwak
* Matrix based on template class, template function
* which means, this library contains dynamic allocation free matrix
*/
#ifndef _MY_DENSE_MATRIX_H_
#define _MY_DENSE_MATRIX_H_

#include <iostream>
#include <iomanip>

// Dense Matrix class for Numerical Analysis
// using template to create matrix in generic programming

template<typename T, int Row, int Colum>
class Dense_Matrix
{

public:
	T data_table[Row][Colum]; // number of row and colum is passed by template parameter,
							  // don't have to allocate memory dynamically

	// constructor & destructor
	Dense_Matrix();        // default constructor, no initializing
	Dense_Matrix(T scala); // matrix with all contents are initialized in scala value, useful for scala multiplication
	Dense_Matrix(const T** input_data); 
	Dense_Matrix(const T input_data[Row][Colum]);
	// can we design this constructor by parameter pack?
	Dense_Matrix(const Dense_Matrix<T, Row, Colum>& mat); // copy constructor

	virtual ~Dense_Matrix(); // don't have to delete memory, 
							 // there  is no memory allocated dynamically

	// io 
	void Set(int row, int col, T value); // set value at data_table at row, colum
	T get(int row, int col); // get value from data_table at row, colum
	void show(); // print the conten of the matrix

	// operator overloading
	// only available for same row and colum
	Dense_Matrix<T, Row, Colum>& operator=(const Dense_Matrix<T, Row, Colum>& mat);
	Dense_Matrix<T, Row, Colum>& operator+=(const Dense_Matrix<T, Row, Colum>& mat);
	Dense_Matrix<T, Row, Colum>& operator-=(const Dense_Matrix<T, Row, Colum>& mat);

	// add & subtract
	// only available for smae row and colum
	friend Dense_Matrix<T, Row, Colum> operator+(const Dense_Matrix<T, Row, Colum>& mat1, const Dense_Matrix<T, Row, Colum>& mat2) {
		T temp[Row][Colum];
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Colum; j++) {
				temp[i][j] = mat1.data_table[i][j] + mat1.data_table[i][j];
			}
		}
		return Dense_Matrix<T, Row, Colum>(temp);
	}
	friend Dense_Matrix<T, Row, Colum> operator-(const Dense_Matrix<T, Row, Colum>& mat1, const Dense_Matrix<T, Row, Colum>& mat2) {
		T temp[Row][Colum];
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Colum; j++) {
				temp[i][j] = mat1.data_table[i][j] - mat1.data_table[i][j];
			}
		}
		return Dense_Matrix<T, Row, Colum>(temp);
	}
	
	// matrix divided by scala value
	template<typename scala_type>
	friend Dense_Matrix<T, Row, Colum> operator/(const Dense_Matrix<T, Row, Colum>& mat, scala_type scala) {
		T temp[Row][Colum];
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Colum; j++) {
				temp[i][j] = mat.data_table[i][j] / scala;
			}
		}
		return Dense_Matrix<T, Row, Colum>(temp);
	}
	
	// etc
	Dense_Matrix<T, Colum, Row> TransPose(); // return transposed matrix
	void eye(); // turn the matrix into Identity matrix 

};

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::Dense_Matrix() {}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::Dense_Matrix(T scala) {

	for (int i = 0; i < Row; i++) {

		for (int j = 0; j < Colum; j++)
			data_table[i][j] = scala;
	}
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::Dense_Matrix(const T** input_data) {

	for (int i = 0; i < Row; i++) {

		for (int j = 0; j < Colum; j++)
			data_table[i][j] = input_data[i][j];
	}
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::Dense_Matrix(const T input_data[Row][Colum]) {

	for (int i = 0; i < Row; i++) {

		for (int j = 0; j < Colum; j++)
			data_table[i][j] = input_data[i][j];
	}
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::Dense_Matrix(const Dense_Matrix<T, Row, Colum>& mat) { // why??????

	for (int i = 0; i < Row; i++) {

		for (int j = 0; j < Colum; j++)
			data_table[i][j] = mat.data_table[i][j];
	}
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>::~Dense_Matrix() {};

template<typename T, int Row, int Colum>
void Dense_Matrix<T, Row, Colum>::Set(int row, int col, T value) {
	data_table[row - 1][col - 1] = value;
}

template<typename T, int Row, int Colum>
T Dense_Matrix<T, Row, Colum>::get(int row, int col) {
	return data_table[row - 1][col - 1];
}

template<typename T, int Row, int Colum>
void Dense_Matrix<T, Row, Colum>::show() {
	for (int i = 0; i < 3 * Colum; i++)	std::cout << "--";
	std::cout << std::endl;
	std::cout << Row << " X " << Colum << " Matrix." << std::endl;
	std::cout << std::endl;
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Colum; j++) {
			std::cout << std::setw(10) << data_table[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < 3 * Colum; i++)	std::cout << "--";
	std::cout << std::endl;
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>& Dense_Matrix<T, Row, Colum>::operator=(const Dense_Matrix<T, Row, Colum>& mat) {
	for (int i = 0; i < Row; i++) {

		for (int j = 0; j < Colum; j++)
			data_table[i][j] = mat.data_table[i][j];
	}
	return *this;
}
template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>& Dense_Matrix<T, Row, Colum>::operator+=(const Dense_Matrix<T, Row, Colum>& mat) {
	(*this) = (*this) + mat;
	return *this;
}
template<typename T, int Row, int Colum>
Dense_Matrix<T, Row, Colum>& Dense_Matrix<T, Row, Colum>::operator-=(const Dense_Matrix<T, Row, Colum>& mat) {
	(*this) = (*this) - mat;
	return *this;
}

template<typename T, int Row, int Colum>
Dense_Matrix<T, Colum, Row> Dense_Matrix<T, Row, Colum>::TransPose() {
	T temp[Colum][Row];
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Colum; j++) {
			temp[j][i] = data_table[i][j];
		}
	}
	return Dense_Matrix<T, Colum, Row> (temp);
}

template<typename T, int Row, int Colum>
void Dense_Matrix<T, Row, Colum>::eye() {
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Colum; j++) {
			if (i == j) {
				data_table[i][j] = 1;
			}
			else {
				data_table[i][j] = 0;
			}
		}
	}
}


// ----------------------------- Not a member function -------------------------- //
// vector calculation
// by defining this function thorough template, we can calculate multiplication between various type of matrix 
template<typename T, int size1, int size2, int size3>
Dense_Matrix<T, size1, size3> operator*(const Dense_Matrix<T, size1, size2>& mat1, const Dense_Matrix<T, size2, size3>& mat2) {
	T temp_matrix[size1][size3];
	for (int i = 0; i < size1; i++) {
		for (int j = 0; j < size3; j++) {
			T temp_value = 0;
			for (int k = 0; k < size2; k++) {
				temp_value += mat1.data_table[i][k] * mat2.data_table[k][j];
			}
			temp_matrix[i][j] = temp_value;
		}
	}
	return Dense_Matrix<T, size1, size3>(temp_matrix);
}

#endif
