#include "PCA.h"

std::vector<double> JacobiMethod(std::vector<std::vector<double>> Matrix){
    int n = (int)Matrix.size();
    double sum;
    std::vector<double> lambdas(n, 0);
    do{
        double maximum = Matrix[0][1];
        int p = 0,q = 1;
        //寻找非对角线按模最大元素
        for(int i = 0; i < n; ++i){
            for(int j = i+1; j < n; ++j)
                if (maximum < fabs(Matrix[i][j])) {
                    maximum = fabs(Matrix[i][j]);
                    p = i, q = j;
                }
        }
        //确定旋转角度
        double s = (Matrix[q][q] - Matrix[p][p]) / (2 * Matrix[p][q]);
        double t1, t2, t, c, d;
        if (s == 0)
            t = 0;
        else{
            t1 = -s - sqrt(s * s + 1);
            t2 = -s + sqrt(s * s + 1);
            t = (fabs(t1) > fabs(t2)) ? t2 : t1;
        }
        c = 1 / sqrt(1 + t * t);
        d = t / sqrt(1 + t * t);
        //计算Q^TAQ相应元素，并存放在A中
        for(int i = 0; i < n; ++i){
            if(i == p || i == q) continue;
            Matrix[i][p] = Matrix[p][i] = c * Matrix[p][i] - d * Matrix[p][i];
            Matrix[i][q] = Matrix[q][i] = c * Matrix[q][i] - d * Matrix[q][i];
        }
        Matrix[p][p] -= t * Matrix[p][q];
        Matrix[q][q] += t * Matrix[p][q];
        Matrix[p][q] = Matrix[q][p] = 0;
    } while((sum = nonDiagonalElements(Matrix)) > 1e-6);

    for(int i = 0; i < n; ++i)
        lambdas[i] = Matrix[i][i];
    return lambdas;
}


void ComputeU(std::vector<std::vector<double>> A, std::vector<double> &lams,
              std::vector<std::vector<double>> &Matrix){
    //计算AA^T及其特征值,对特征值进行排序
    std::vector<std::vector<double>> AAT = MatrixProduct(A, MatrixTranspose(A));
    lams = JacobiMethod(AAT);
    std::sort(lams.begin(), lams.end(), std::greater<>());
    //TODO
}

double nonDiagonalElements(std::vector<std::vector<double>> Matrix){
    //返回矩阵所有非对角元素的平方和
    int n = (int)Matrix.size();
    double sum = 0;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(i == j) continue;
            sum += Matrix[i][j] * Matrix[i][j];
        }
    }
    return sum;
}

std::vector<std::vector<double>> MatrixTranspose(std::vector<std::vector<double>> Matrix){
    //矩阵转置
    int M_rows = (int)Matrix.size(), M_cols = (int)Matrix[0].size();
    std::vector<std::vector<double>> M_Trans(M_cols, std::vector<double>(M_rows, 0));
    for(int i = 0; i < M_rows; ++i)
        for(int j = 0; j <M_cols; ++j)
            M_Trans[j][i] = Matrix[i][j];
    return M_Trans;
}

std::vector<std::vector<double>> MatrixProduct(std::vector<std::vector<double>> A,
                                               std::vector<std::vector<double>> B){
    int A_rows = (int)A.size(), A_cols = (int)A[0].size();
    int B_rows = (int)B.size(), B_cols = (int)B[0].size();
    //先判断是否能相乘
    if(A_cols != B_rows){
        std::cout << "Matrices cannot be multiplied!"<<std::endl;
        exit(255);
    }
    //根据定义相乘
    std::vector<std::vector<double>> Product(A_rows, std::vector<double>(B_cols, 0));
    for(int i = 0;i < A_rows; ++i){
        for(int j = 0; j < B_cols; ++j){
            double sum = 0;
            for(int k = 0; k < A_cols; ++k) sum += A[i][k] * B[k][j];
            Product[i][j] = sum;
        }
    }
    return Product;
}