#include <iostream>
#include "LinearEquations.h"

int main(){
    //接收参数
    double epsilon, a, h;
    int n;
    std::cout << "请输入epsilon,a,n: ";
    std::cin >> epsilon >> a >> n;
    h = 1.0/n;
    //创建系数矩阵和常数项矩阵
    auto coeff = InitMatrix(n-1, epsilon, h);
    std::vector<double> b(n-1,a*h*h);
    b[n-2] -= epsilon+h;
    //两种方法分别计算
    auto solve1 = GaussElimination(n-1, coeff, b);
    auto solve2 = Gauss_Seidel(n-1, coeff, b, 1e-15);
    auto Tvalue = ComputeExact(n,epsilon,a);
    //打印结果
    std::cout << "Gauss列主元结果: ";
    for(auto item: solve1) std::cout << item << " " ;
    std::cout << "\nGauss_Seidel结果: ";
    for(auto item: solve2) std::cout << item << " " ;
    std::cout << "\n真实值: ";
    for(auto item: Tvalue) std::cout << item << " " ;

    return 0;
}