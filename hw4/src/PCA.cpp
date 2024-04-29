#include "PCA.h"

void JacobiMethod(std::vector<std::vector<double>> Matrix, std::vector<double> &lambdas, int debug){
    int n = (int)Matrix.size(), times = 0;
    double sum = nonDiagonalElements(Matrix);
    while(sum > 1e-6){
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
            t = 1;
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
            double M_pi = Matrix[p][i], M_qi = Matrix[q][i];
            Matrix[i][p] = Matrix[p][i] = c * M_pi - d * M_qi;
            Matrix[i][q] = Matrix[q][i] = c * M_qi + d * M_pi;
        }
        Matrix[p][p] -= t * Matrix[p][q];
        Matrix[q][q] += t * Matrix[p][q];
        Matrix[p][q] = Matrix[q][p] = 0;
        if(debug)
            std::cout << "第"<< times <<  "次迭代， 非对角元素平方和： " << sum << std::endl;
        times++;
        sum = nonDiagonalElements(Matrix);
    }
    std::cout << "最终非对角元素平方和： " << sum << std::endl;
    for(int i = 0; i < n; ++i)
        lambdas[i] = Matrix[i][i];
}


void ComputeU(std::vector<std::vector<double>> A, std::vector<double> &lams,
              std::vector<std::vector<double>> &U){
    int n = (int)A.size();
    //计算AA^T及其特征值,对特征值进行排序
    std::vector<std::vector<double>> AAT = MatrixProduct(A, MatrixTranspose(A));
    JacobiMethod(AAT, lams,1);
    std::sort(lams.begin(), lams.end(), std::greater<>());
    //求矩阵U
    std::vector<double> fea(n, 0);
    std::vector<double> zero(n, 0);
    std::vector<std::vector<double>> AAT_copy(n, std::vector<double>(n,0));
    for(int i = 0; i < n; ++i){
        AAT_copy = AAT;
        for(int j = 0; j < n; ++j)
            AAT_copy[j][j] -= lams[i];
        fea = GaussElimination(AAT_copy, zero);
        for(int j = 0; j < n; ++j)
            U[j][i] = fea[j];
    }
}

double nonDiagonalElements(std::vector<std::vector<double>> Matrix){
    //返回矩阵所有非对角元素的平方和
    int n = (int)Matrix.size();
    double sum = 0;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(i == j)
                continue;
            else
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


std::vector<double> GaussElimination(std::vector<std::vector<double>> coeff,
                                     std::vector<double> b){
    //高斯消元法求解方程组
    int n = (int)coeff.size();
    std::vector<double> y(n, 0);
    //列主元消元
    for(int i = 0; i < n-1; ++i){
        //寻找并交换列主元
        int k=i;
        for(int j = i+1; j < n; ++j)
            if(fabs(coeff[k][i]) < fabs(coeff[j][i]))
                k = j;
        std::vector<double> temp = coeff[i];
        coeff[i] = coeff[k], coeff[k] = temp;
        //消元
        for(int j = i+1; j < n; ++j){
            double t = coeff[j][i] / coeff[i][i];
            for(int h=i; h < coeff[j].size(); ++h)
                coeff[j][h] -= t * coeff[i][h];
        }
    }
    //std::cout << "处理前： ", putVector(y);
    int num = 0;
    for(int i = 0; i < n; ++i){
        if(InfiniteModulus(coeff[i]) < 1e-3) num++;
    }

    for(int i = 0; i < num; ++i)
        y[n-1-i] = 1;
    //std::cout << "处理num后： ", putVector(y);
    //std::cout<< "num: " <<num<<std::endl;
    for(int i = 0; i < n-num; ++i)
        for(int j = 0; j < num; ++j)
            b[i] -= coeff[i][n-1-j] * y[n-1-j];
    //回代
    for(int i=n-num-1; i>=0; --i){
        y[i] = b[i] / coeff[i][i];
        for(int j=i-1; j >= 0; --j)
            b[j] -= coeff[j][i] * y[i];
    }
    //归一化
    double sum = 0;

    for(auto item: y){
        sum += item * item;
    }
    std::cout << std::endl;
    for(int i = 0; i < n; ++i)
        y[i] /= sqrt(sum);
    return y;
}

std::vector<double> LU(std::vector<double> b, std::vector<std::vector<double>> A){
    int n = (int)A[0].size(), m = (int)A.size();
    for(int i=0;i<m-n;i++)  A.pop_back();
    std::vector<double> y(n, 0);
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0));
    std::vector<double> X_m(n,1);
    for (int i = 0; i < n; i++)
        L[i][i] = 1;
    for (int k = 0; k < n; k++){
        for (int j = k; j < n; j++){
            U[k][j] = A[k][j];
            for (int r = 0; r < k; r++)
                U[k][j] -= L[k][r] * U[r][j];
        }
        for (int i = k + 1; i < n; i++){
            L[i][k] = A[i][k];
            for (int r = 0; r < k; r++)
                L[i][k] -= L[i][r] * U[r][k];
            L[i][k] /= U[k][k];
        }
    }
    for (int i = 0; i < n; i++){
        X_m[i] = b[i];
        for (int j = 0; j < i; j++)
            X_m[i] -= L[i][j] * X_m[j];
    }
    for (int i = n - 1; i >= 0; i--){
        y[i] = X_m[i];
        for (int j = i + 1; j < n; j++)
            y[i] -= U[i][j] * y[j];
        y[i] /= U[i][i];
    }
    return y;
}

double determinant(std::vector<std::vector<double>> A){
    //计算方阵A的行列式值
    int n = (int)A.size(), m = (int)A[0].size();
    //先判断是不是方阵
    if(n != m){
        std::cout << "Do u know the story of Jiafeng Tang?(fog)"<<std::endl;
        exit(255);
    }
    //用高斯消元法，最后将主对角线元素相乘
    for(int i=0; i<n; ++i){
        //寻找并交换列主元
        int k=i;
        for(int j=i+1; j<n; ++j)
            if(fabs(A[k][i]) < fabs(A[k][j]))
                k = j;
        std::vector<double> temp = A[i];
        A[i] = A[k], A[k] = temp;
        //消元
        for(int j=i+1; j<n; ++j){
            double t = A[j][i] / A[i][i];
            for(int h=i; h < A[j].size(); ++h)
                A[j][h] -= t * A[i][h];
        }
    }
    double product = 1;
    for(int i = 0; i < n; ++i)
        product *= A[i][i];
    return product;
}
//字符串分割
std::vector<std::string> split(std::string s, const std::string &delimiter) {
    std::vector<std::string> res;
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        res.push_back(token);
        s = s.substr(pos + delimiter.length());
    }
    res.push_back(s);
    return res;
}
//计算两向量内积
double vectorProduct(std::vector<double> a, std::vector<double> b){
    if(a.size() != b.size()){
        std::cout << "vectors have different divs!" << std::endl;
        exit(255);
    }
    double sum = 0;
    for(int i = 0; i < a.size(); ++i)
        sum += a[i] * b[i];
    return sum;
}
//在控制台打印向量和矩阵
void putVector(std::vector<double> v){
    for(auto item: v)
        std::cout << item << ", ";
    std::cout << std::endl;
}
void putMatrix(std::vector<std::vector<double>> matrix){
    for(auto item: matrix)
        putVector(item);
}

double InfiniteModulus(std::vector<double> v){
    double max = fabs(v[0]);
    for(auto item: v)
        if (max < fabs(item)) max = fabs(item);
    return max;
}