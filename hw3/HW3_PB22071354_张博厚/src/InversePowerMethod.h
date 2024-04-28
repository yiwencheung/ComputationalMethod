#include <iostream>
#include <vector>
#include <cmath>

//LU分解：输入系数矩阵A和L，U的引用
void Doolittle(std::vector<std::vector<double>>,
               std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
//回代求解
std::vector<double> LU_Solve(std::vector<std::vector<double>>,
                             std::vector<std::vector<double>>, std::vector<double>);
//反幂法迭代
void InversePowerMethod(std::vector<std::vector<double>> A);

//估计特征值：寻找绝对值最大的分量
double find_max(std::vector<double> x);




