#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>

#define FILE_PATH "D:\\Projects\\ComputationalMethod\\hw4\\res\\iris.txt"
#define SAVE_PATH0 "D:\\Projects\\ComputationalMethod\\hw4\\res\\label0.txt"
#define SAVE_PATH1 "D:\\Projects\\ComputationalMethod\\hw4\\res\\label1.txt"
#define SAVE_PATH2 "D:\\Projects\\ComputationalMethod\\hw4\\res\\label2.txt"

void ComputeU(std::vector<std::vector<double>>, std::vector<double>&, std::vector<std::vector<double>>&);
//Jacobi算法求矩阵特征值
void JacobiMethod(std::vector<std::vector<double>>, std::vector<double>&, int);
//计算矩阵所有非对角元平方和
double nonDiagonalElements(std::vector<std::vector<double>>);
//计算矩阵转置
std::vector<std::vector<double>> MatrixTranspose(std::vector<std::vector<double>>);
//计算矩阵乘积
std::vector<std::vector<double>> MatrixProduct(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
//高斯消元法
std::vector<double> GaussElimination(std::vector<std::vector<double>>, std::vector<double>);
std::vector<double> LU(std::vector<double>,std::vector<std::vector<double>>);
//计算行列式值
double determinant(std::vector<std::vector<double>>);
//经典字符串分割算法
std::vector<std::string> split(std::string, const std::string &);
//两向量内积
double vectorProduct(std::vector<double>, std::vector<double>);
//在控制台打印向量和矩阵
void putVector(std::vector<double>);
void putMatrix(std::vector<std::vector<double>>);
//计算向量的无穷范数
double InfiniteModulus(std::vector<double>);