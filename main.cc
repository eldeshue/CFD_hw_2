/* 2022.10.03
 * 전산 유체 역학 과제 1
 * y = x^3 함수의 x=1에서의 미분값 3 에 대하여,
 * 1계 정확도의 전방 차분법과 2계 정확도의 중심 차분법에 의해서 발생하는 에러에 대해서 비교한다.
 * compile option: g++ -o main main.cc -std=c++17 -lboost_iostreams -lboost_system -lboost_filesystem
 */

#include<iostream>
#include<iomanip>
#include<functional>
#include<cmath>
#include<vector>
#include<map>
#include"gnuplot-iostream.h" // library for connecting gnuplot

int main(){
	std::cout << std::setprecision(16);
  // set Gnuplot object;
	Gnuplot gp;			
	// set container for the graph
  //std::vector<std::pair<long double,long double>> data;
	
	
	// make a plot 
	// y value is   and x value is 
	// setting option for the gnuplot
	gp << "set xlabel 'x축 이름' font 'Helvetica, 12'\n" << 
    "set ylabel 'y축 이름' font 'Helvetica, 12'" << std::endl;
	gp << "set tics font 'Helvetica, 10'" << std::endl;
	gp << "set title '전체 그래프 이름' font 'Helvetica, 13'" << std::endl;
	gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" << 
    "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
  //gp << "set xrange [-1:30]\nset yrange [-15:2]\n";
	gp << "set autoscale \n"; 
  gp << "set grid \n";
  //gp << "plot" << gp.file1d(sided_arr) << "with points title 'sided',"
  //<< gp.file1d(central_arr) << "with points title 'central'" << std::endl;
  gp << "plot '-' with lines title 'sided', '-' with lines title 'central'\n";
  // send data to show
	//gp.send1d(data);
}
