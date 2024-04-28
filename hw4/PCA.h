#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

void ComputeU(std::vector<std::vector<double>>, std::vector<double>&, std::vector<std::vector<double>>&);
//Jacobi算法求矩阵特征值
std::vector<double> JacobiMethod(std::vector<std::vector<double>> );
//计算矩阵所有非对角元平方和
double nonDiagonalElements(std::vector<std::vector<double>>);
//计算矩阵转置
std::vector<std::vector<double>> MatrixTranspose(std::vector<std::vector<double>>);
//计算矩阵乘积
std::vector<std::vector<double>> MatrixProduct(std::vector<std::vector<double>>, std::vector<std::vector<double>>);