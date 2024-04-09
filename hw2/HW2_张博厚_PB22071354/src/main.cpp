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
    //��������ֵ�뾫ȷֵ�����
    auto loss1 = InfinitePara(solve1, Tvalue);
    auto loss2 = InfinitePara(solve2, Tvalue);
    auto diff = InfinitePara(solve1,solve2);
    //��ӡ���
    std::cout << "Gauss����Ԫ���: ";
    for(auto item: solve1) std::cout << item << " " ;
    std::cout << "\nGauss_Seidel���: ";
    for(auto item: solve2) std::cout << item << " " ;
    std::cout << "\n��ʵֵ: ";
    for(auto item: Tvalue) std::cout << item << " " ;
    std::cout << "\n���������ֵ������:" << diff;
    std::cout << "\nGauss����Ԫ������ʵֵ������:" << loss1 ;
    std::cout << "\nGauss-Seidel����������ʵֵ������:" << loss2 << std::endl;
    return 0;
}