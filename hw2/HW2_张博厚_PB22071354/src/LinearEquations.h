#include <iostream>
#include <vector>
#include <cmath>

std::vector<std::vector<double>> InitMatrix(int n, double epsilon, double h);   //创建系数矩阵

std::vector<double> GaussElimination(int n,                                     //高斯列主元消元法
                    std::vector<std::vector<double>> coeff, std::vector<double> b);

std::vector<double> Gauss_Seidel(int n, std::vector<std::vector<double>> coeff,
                                 std::vector<double> b, double accuracy);       //高斯-赛德尔迭代

std::vector<double> ComputeExact(int n, double epsilon, double a);                     //计算真实值

double InfinitePara(std::vector<double> a, std::vector<double> b);              //计算两向量差的无穷范数



