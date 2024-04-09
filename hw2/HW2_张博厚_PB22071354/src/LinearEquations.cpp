#include "LinearEquations.h"

std::vector<std::vector<double>> InitMatrix(int n, double epsilon, double h){
    //创建矩阵
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));
    //设置系数矩阵
    matrix[0][0] = -(2*epsilon+h), matrix[0][1] = epsilon+h;
    matrix[n-1][n-2] = epsilon, matrix[n-1][n-1] =  -(2*epsilon+h);
    for(int i=1; i<n-1; ++i){
        matrix[i][i-1] = epsilon;
        matrix[i][i] = -(2*epsilon+h);
        matrix[i][i+1] = epsilon+h;
    }
    return matrix;
}

std::vector<double> GaussElimination(int n,
                    std::vector<std::vector<double>> coeff, std::vector<double> b){
    //arg:系数矩阵 b：常数项矩阵
    //列主元消元
    for(int i=0; i<n; ++i){
        //寻找并交换列主元
        int k=i;
        for(int j=i+1; j<n; ++j)
            if(fabs(coeff[k][i]) < fabs(coeff[k][j]))
                k = j;
        std::vector<double> temp = coeff[i];
        coeff[i] = coeff[k], coeff[k] = temp;
        double temp_b = b[i];
        b[i] = b[k], b[k] = temp_b;
        //消元
        for(int j=i+1; j<n; ++j){
            double t = coeff[j][i] / coeff[i][i];
            for(int h=i; h < coeff[j].size(); ++h)
                coeff[j][h] -= t * coeff[i][h];
            b[j] -= t*b[i];
        }
    }
    //回代
    for(int i=n-1; i>=0; --i){
        for(int j=i+1; j<n; ++j)
            b[i] -= coeff[i][j]*b[j];
        b[i] /= coeff[i][i];
    }
    return b;
}

std::vector<double> Gauss_Seidel(const int n, std::vector<std::vector<double>> coeff,
                                 std::vector<double> b,const double accuracy){
    //设置迭代初始值为0
    std::vector<double> x1(n,0);
    std::vector<double> x2(n,1);
    while(InfinitePara(x1,x2) > accuracy){
        x1 = x2;
        //由x(k)计算x(k+1)
        for(int i=0; i<n; ++i){
            double sum=0;
            for(int j=0; j<n; ++j)
                sum += coeff[i][j]*x2[j];
            x2[i] = (b[i] - sum + coeff[i][i]*x2[i])/coeff[i][i];
        }
    }
    return x2;
}

std::vector<double> ComputeExact(int n, double epsilon, double a){
    std::vector<double> TValue(n-1,0);
    double h = 1.0/n, x;
    for(int i=1; i<n; ++i){
        x = i*h;
        TValue[i-1] = a*x + (1-a)*(1 - exp(-x/epsilon))/(1 - exp(-1/epsilon));
    }
    return TValue;
}

double InfinitePara(std::vector<double> a, std::vector<double> b){
    if(a.size() != b.size()){
        std::cout << "两向量维度不同!" << std::endl;
        exit(0);
    }
    double max = fabs(a[0]-b[0]);
    for(int i=1; i < a.size(); ++i)
        if(max < fabs(a[i]-b[i]))
            max = fabs(a[i]-b[i]);
    return max;
}