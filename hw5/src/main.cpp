#include "Spline.h"

int main(){
    /*计算三次样条插值
     * 使用大M法，取自然边界条件。
    */
    //所用到的变量
    int h = 1;      //步长
    int i = 0;
    double lambda = 0.5, mu = 0.5;
    std::fstream fin;
    std::string rdline;
    std::vector<std::vector<double>> points(POINTNUM, std::vector<double>(DIV, 0));
    std::vector<double> M(POINTNUM, 0);
    std::vector<std::vector<double>> S(POINTNUM - 1, std::vector<double>(POWER + 1, 0));
    //读取并分割数据
    fin.open(FILE_PATH, std::ios::in);
    while(std::getline(fin, rdline)){
        std::vector<std::string> data = split(rdline, "\t");
        points[i][0] = std::stod(data[0]);
        points[i][1] = std::stod(data[1]);
        ++i;
    }
    fin.close();
    //求解M
    M = CubicSpline_natural(lambda, mu, h, points);
    //由M求解多项式系数
    S = SolveS(M, points, h);
    //输出多项式
    PolyDisp(S, points);

    //更换点的位置重新计算
    std::vector<double> M_plus(POINTNUM, 0);
    std::vector<std::vector<double>> S_plus(POINTNUM - 1, std::vector<double>(POWER + 1, 0));
    points[9][1] = 10;
    M_plus = CubicSpline_natural(lambda, mu, h, points);
    S_plus = SolveS(M_plus, points, h);
    std::cout << "\n更改第十个点后，结果为:" << std::endl;
    PolyDisp(S_plus, points);

    //比较前后两次结果
    std::cout << "更改前后，多项式的系数变化为：\n";
    for(int j = 0; j < S.size(); ++j){
        std::cout << S[j][0] - S_plus[j][0] << ", ";
        std::cout << S[j][1] - S_plus[j][1] << ", ";
        std::cout << S[j][2] - S_plus[j][2] << ", ";
        std::cout << S[j][3] - S_plus[j][3] << std::endl;
    }

    //保存结果
    std::fstream fout0, fout1;
    fout0.open(SAVE_PATH_0, std::ios::out);
    fout1.open(SAVE_PATH_1, std::ios::out);
    for(auto s: S)
        fout0 << s[0] << " " << s[1] << " " << s[2] << " " << s[3] <<std::endl;
    for(auto s: S_plus)
        fout1 << s[0] << " " << s[1] << " " << s[2] << " " << s[3] <<std::endl;
    fout0.close();
    fout1.close();

    return 0;
}