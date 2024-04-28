#include "PCA.h"


int main(){
    //创建随机矩阵
    int A_rows = 4, A_cols = 3;
    std::vector<std::vector<double>> A(A_rows, std::vector<double>(A_cols, 0));
    std::default_random_engine rand_e;
    std::uniform_real_distribution<double> rand_uni(0, 1.0);

    rand_e.seed(114514);
    for(int i = 0; i < A_rows; ++i){
        for(int j = 0; j < A_cols; ++j)
           A[i][j] = rand_uni(rand_e);
    }
    //计算AA^T的特征值和特征向量
    std::vector<double> U_lams(A_rows,0);
    std::vector<std::vector<double>> U(A_rows,std::vector<double> (A_rows,0));



}