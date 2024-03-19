#include <iostream>
#include "LinearEquations.h"

int main(){
    //���ղ���
    double epsilon, a, h;
    int n;
    std::cout << "������epsilon,a,n: ";
    std::cin >> epsilon >> a >> n;
    h = 1.0/n;
    //����ϵ������ͳ��������
    auto coeff = InitMatrix(n-1, epsilon, h);
    std::vector<double> b(n-1,a*h*h);
    b[n-2] -= epsilon+h;
    //���ַ����ֱ����
    auto solve1 = GaussElimination(n-1, coeff, b);
    auto solve2 = Gauss_Seidel(n-1, coeff, b, 1e-15);
    auto Tvalue = ComputeExact(n,epsilon,a);
    //��ӡ���
    std::cout << "Gauss����Ԫ���: ";
    for(auto item: solve1) std::cout << item << " " ;
    std::cout << "\nGauss_Seidel���: ";
    for(auto item: solve2) std::cout << item << " " ;
    std::cout << "\n��ʵֵ: ";
    for(auto item: Tvalue) std::cout << item << " " ;

    return 0;
}