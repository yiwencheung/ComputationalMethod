#include "PCA.h"


int main(){
    //part 1
    std::cout << "----------PART I-----------" << std::endl;
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
    std::cout << "生成的随机矩阵为:\n";
    putMatrix(A);
    //计算AA^T的特征值和特征向量
    std::vector<double> U_lams(A_rows,0);
    std::vector<std::vector<double>> U(A_rows,std::vector<double> (A_rows,0));
    ComputeU(A, U_lams, U);
    std::cout << "AA^T的特征值为：";
    putVector(U_lams);
    std::cout << "矩阵U为:\n";
    putMatrix(U);
    //计算A^TA的特征值并排序
    std::vector<std::vector<double>> ATA = MatrixProduct(MatrixTranspose(A), A);
    std::vector<double> V_lams(A_cols, 0);
    JacobiMethod(ATA, V_lams, 0);
    std::sort(V_lams.begin(), V_lams.end(), std::greater<double>());
    //计算矩阵V
    std::vector<std::vector<double>> Singular(A_rows, std::vector<double> (A_cols,0));
    for(int i=0;i < A_rows; i++)    Singular[i][i]=sqrt(V_lams[i]);
    std::cout << "奇异值矩阵为:\n";
    putMatrix(Singular);
    std::vector<std::vector<double>> temp = MatrixTranspose(MatrixProduct(U, Singular));
    std::vector<std::vector<double>> V(A_cols, std::vector<double> (A_cols,0));
    std::vector<double> Y(A_cols, 0);
    for(int i = 0; i < A_cols; ++i){
        Y = LU(temp[i], A);
        for(int j = 0; j < A_cols; ++j)
            V[j][i] = Y[j];
    }
    std::cout << "矩阵V为:\n";
    putMatrix(V);
    //结果分析
    std::vector<std::vector<double>> AAT = MatrixProduct(A, MatrixTranspose(A));
    std::vector<std::vector<double>> Ana(A_rows, std::vector<double>(A_rows, 0));
    for(int i = 0; i < A_rows; ++i){
        Ana = AAT;
        for(int j = 0; j < A_rows; ++j)
            Ana[j][j] -= U_lams[i];
        std::cout << "det(AAT - lambda" <<i <<" I) = " << determinant(Ana) << std::endl;
    }

    //part 2
    std::cout << "----------PART II-----------" << std::endl;
    std::ifstream fs(FILE_PATH, std::ios::in);
    if(!fs.is_open()){
        std::cout << "open file failed!" << std::endl;
        exit(255);
    }
    //将数据存入矩阵,列向量形式
    int data_num = 150, data_div = 4;
    std::vector<std::vector<double>> B(data_div, std::vector<double>(data_num,0));
    std::string rdline;
    int label[3] = {0}, j = 0;
    while(std::getline(fs, rdline)){
        std::vector<std::string> datas = split(rdline, ",");
        for(int i = 0; i < 5; ++i){
            if(i != 4)
                B[i][j] = std::stod(datas[i]);
            else
                label[std::stoi(datas[i])]++;
        }
        ++j;
    }
    //数据去中心化
    std::vector<double> average(data_div, 0);
    for(int i = 0; i < data_div; ++i){
        average[i] = std::accumulate(B[i].begin(), B[i].end(), 0.00) / data_num;
        for(int j = 0; j < data_num; ++j)
            B[i][j] -= average[i];
    }
    //计算协方差矩阵1/m BB^T
    std::vector<std::vector<double>> BBT = MatrixProduct(B, MatrixTranspose(B));
    for(int i = 0; i < data_div; ++i)
        for(int j = 0; j < data_div; ++j)
            BBT[i][j] /= data_num;
    std::cout << "协方差矩阵为:\n";
    putMatrix(BBT);
    //求协方差矩阵的特征值
    std::vector<double> B_lams(data_div, 0);
    JacobiMethod(BBT, B_lams, 1);
    std::sort(B_lams.begin(), B_lams.end(), std::greater<double>());
    //求解基向量
    std::vector<double> e1(data_div, 0), e2(data_div, 0), zero(data_div, 0);
    std::vector<std::vector<double>> BBT_copy(data_div, std::vector<double>(data_div, 0));
    for(int i = 0; i < 2; ++i){
        BBT_copy = BBT;
        for(int j = 0; j < data_div; ++j)
            BBT_copy[j][j] -= B_lams[i];
        if(i == 0)
            e1 = GaussElimination(BBT_copy, zero);
        else
            e2 = GaussElimination(BBT_copy, zero);
    }
    //保存数据
    double label0_points[label[0]][2], label1_points[label[1]][2], label2_points[label[2]][2];
    std::vector<std::vector<double>> BT = MatrixTranspose(B);
    std::fstream l0, l1, l2;
    l0.open(SAVE_PATH0, std::ios::out);
    l1.open(SAVE_PATH1, std::ios::out);
    l2.open(SAVE_PATH2, std::ios::out);
    for(int i = 0; i < label[0]; ++i)
        l0 << vectorProduct(BT[i], e1) << "," << vectorProduct(BT[i], e2) << std::endl;
    for(int i = label[0]; i < label[0] + label[1]; ++i)
        l1 << vectorProduct(BT[i], e1) << "," << vectorProduct(BT[i], e2) << std::endl;
    for(int i = label[0] + label[1]; i < label[0] + label[1] + label[2]; ++i)
        l2 << vectorProduct(BT[i], e1) << "," << vectorProduct(BT[i], e2) << std::endl;
    l0.close();
    l1.close();
    l2.close();

    return 0;
}
