#include "InversePowerMethod.h"

double find_max(std::vector<double> x) {
    double max = fabs(x[0]);
    for (int i = 1; i < x.size(); ++i)
        if (fabs(x[i]) > max)
            max = fabs(x[i]);
    return max;
}

void Doolittle(std::vector<std::vector<double>> A,
               std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U){
    int n = A.size();
    for(int i = 0; i < n; ++i){
        //����U�ĵ�i��Ԫ��
        for (int j = i; j < n; ++j){
            double sum = 0;
            for(int r = 0; r < i; ++r)
                sum += L[i][r]*U[r][j];
            U[i][j] = A[i][j] - sum;
        }
        //����L�ĵ�i��Ԫ��
        for (int j = i+1; j < n; ++j) {
            double sum = 0;
            for(int r = 0; r < i; ++r)
                sum += L[j][r]*U[r][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
}

std::vector<double> LU_Solve(std::vector<std::vector<double>> L,
                             std::vector<std::vector<double>> U, std::vector<double> b){
    int n = b.size();
    std::vector<double> y(n,0), x(n,0);
    //�ⷽ����LY = b;
    for(int i = 0; i < n; ++i){
        double sum = 0;
        for(int j = 0; j < i; ++j)
            sum += L[i][j]*y[j];
        y[i] = b[i] - sum;
    }
    //�ⷽ����UX = Y
    for(int i = n-1; i >= 0; --i){
        double sum = 0;
        for(int j = i+1; j < n; ++j)
            sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / (double)U[i][i];
    }
    return x;
}

void InversePowerMethod(std::vector<std::vector<double>> A){
    int n = A.size();
    //��ʼ��L��U��X��Y
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0)), U(n, std::vector<double>(n, 0));
    std::vector<double> x(n,1), y(n,0);
    for(int i = 0; i < n; ++i)          //��ʼ��LΪ�Խ��� UΪ�����
        for(int j = 0; j < n; ++j){
            L[i][j] = (i==j) ? 1 : 0;
            U[i][j] = 0;
        }
    //��LU�ֽ�
    Doolittle(A,L,U);
    //����
    double lambda0 = 0, lambda1 = 1;
    int times = 0;      //ѭ������
    while(fabs(lambda0 - lambda1) >= 1e-5){
        lambda0 = lambda1;
        for(int i = 0; i < y.size(); i++)
            y[i] = x[i] / lambda0;
        //�������ѭ����ֵ
        std::cout << "X" << times << ": ";
        for(auto xi: x)
            std::cout << xi << ", ";
        std::cout << "\nY" << times << ": ";
        for(auto yi: y)
            std::cout << yi << ", ";
        std::cout << "\n ��������ֵ��" << lambda0 << std::endl;
        //����X
        x = LU_Solve(L, U, y);
        lambda1 = find_max(x);
        times ++;
    }
}
