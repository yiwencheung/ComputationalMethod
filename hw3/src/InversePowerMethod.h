#include <iostream>
#include <vector>
#include <cmath>

//需要处理double和int两种数据类型，使用功能模版
//template <typename T>   T max(T a, T b);
//LU分解：输入系数矩阵A和L，U的引用
//template <typename T, typename B>   void Doolittle(T A, T &L, T &U);
//回代求解
//template <typename T, typename N, typename B>   N LU_Solve(T L, T U, B x);
//反幂法迭代
//template <typename T, typename N, typename B>   void InversePowerMethod(T A);

//估计特征值：寻找绝对值最大的分量
//template <typename N, typename B> B find_max(N x);

template <typename T>
T abs(T a) {
    return a>0 ? a : -a;
}

template <typename N, typename B>
B find_max(N x) {
    B max = abs<B>(x[0]);
    for (int i = 1; i < x.size(); ++i)
        if (abs<B>(x[i]) > max)
            max = abs<B>(x[i]);
    return max;
}

template <typename T, typename B>       //T为矩阵类型，B为矩阵中元素类型
void Doolittle(T A, T &L, T &U){
    int n = A.size();
    for(int i = 0; i < n; ++i){
        //计算U的第i行元素
        for (int j = i; j < n; ++j){
            B sum = 0;
            for(int r = 0; r < i; ++r)
                sum += L[i][r]*U[r][j];
            U[i][j] = A[i][j] - sum;
        }
        //计算L的第i列元素
        for (int j = i+1; j < n; ++j) {
            B sum = 0;
            for(int r = 0; r < i; ++r)
                sum += L[j][r]*U[r][i];
            L[j][i] = (A[j][i] - sum) / (double)U[i][i];      //本题中有double和int两种类型，防止int类型做除法出现精度损失
        }
    }
}

template <typename T, typename N, typename B>
N LU_Solve(T L, T U, N b){
    N y = b, x = b;    //确保size相同
    int n = b.size();
    //解方程组LY = b;
    for(int i = 0; i < n; ++i){
        B sum = 0;
        for(int j = 0; j < i; ++j)
            sum += L[i][j]*y[j];
        y[i] = b[i] - sum;
    }
    //解方程组UX = Y
    for(int i = n-1; i >= 0; --i){
        B sum = 0;
        for(int j = i+1; j < n; ++j)
            sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / (double)U[i][i];
    }
    return x;
}

template <typename T, typename N, typename B>
void InversePowerMethod(T A){
    int n = A.size();
    //初始化L，U，X，Y
    T L = A, U = A;
    N x,y;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j){
            L[i][j] = (i==j) ? 1 : 0;
            U[i][j] = 0;
        }
    x.resize(n,1);
    //做LU分解
    Doolittle<T,B>(A,L,U);
    //迭代
    B lambda0 = 0, lambda1 = 1;
    while(abs<B>(lambda0 - lambda1) >= 1e-5){
        lambda0 = lambda1;
        for(int i = 0; i < y.size(); i++)
            y[i] = x[i] / lambda0;
        x = LU_Solve<T, N, B>(L, U, y);
        lambda1 = find_max<N, B>(x);
    }
}



