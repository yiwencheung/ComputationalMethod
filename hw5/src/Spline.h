#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#define POINTNUM    21
#define DIV         2
#define POWER       3
#define FILE_PATH   "D:\\Projects\\ComputationalMethod\\hw5\\res\\point.txt"
#define SAVE_PATH_0 "D:\\Projects\\ComputationalMethod\\hw5\\res\\para0.txt"
#define SAVE_PATH_1 "D:\\Projects\\ComputationalMethod\\hw5\\res\\para1.txt"

std::vector<std::string> split(std::string, const std::string &);
std::vector<double> CubicSpline_natural(double, double, int, std::vector<std::vector<double>>);
double InfinitePara(std::vector<double>, std::vector<double>);
std::vector<double> Gauss_Seidel(int, std::vector<std::vector<double>>,
                                 std::vector<double>,double);
std::vector<std::vector<double>> SolveS(std::vector<double>, std::vector<std::vector<double>>, int);
void PolyDisp(std::vector<std::vector<double>>, std::vector<std::vector<double>>);