#include "InversePowerMethod.h"

int main(){
    //初始化系数矩阵
    std::vector<std::vector<double>>    A2  = {{4, -1 , 1, 3}, {16, -2, -2, 5}, {16, -3, -1, 7}, {6, -4, 2, 9}};
    std::vector<std::vector<double>> A1 (5, std::vector<double>(5,0));
    for(int i = 9, k = 0; i > 4; --i, ++k)
        for (int j = 0; j < 5; ++j)
            A1[k][j] = 1.0 / (i-j);
    //分别求解
    InversePowerMethod(A1);
    InversePowerMethod(A2);
    return 0;
}