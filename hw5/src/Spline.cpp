#include "Spline.h"

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

std::vector<double> CubicSpline_natural(double lambda, double mu,int h,  std::vector<std::vector<double>> points){
    int n = (int)points.size() - 1;
    std::vector<double> d(n-1, 0);
    std::vector<double> M(n-1, 0);
    std::vector<std::vector<double>> coeff(n-1, std::vector<double>(n-1,0));
    //设置系数矩阵
    for(int i = 0; i < n-1; ++i){
        d[i] = 6 * (points[i+2][1] - 2 * points[i+1][1] + points[i][1]) / (2 * h * h);
        coeff[i][i] = 2;
        if(i != 0)
            coeff[i][i-1] = mu;
        if(i != n-2)
            coeff[i][i+1] = lambda;
    }
    //用高斯-赛德尔迭代求解
    M = Gauss_Seidel(n-1, coeff, d, 1e-6);
    M.insert(M.begin(), 0);
    M.push_back(0);
    return M;
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

std::vector<std::vector<double>> SolveS(std::vector<double> M,
                                        std::vector<std::vector<double>> Points, int h){
    std::vector<std::vector<double>> S(POINTNUM - 1, std::vector<double>(POWER + 1, 0));
    for(int i = 0; i < POINTNUM-1; ++i){
        S[i][0] = (M[i+1] - M[i]) / (6.0 * h);
        S[i][1] = (M[i] * Points[i+1][0] - M[i+1] * Points[i][0]) / (2.0 * h);
        S[i][2] = ((M[i+1] * pow(Points[i][0], 2) - M[i] * pow(Points[i+1][0], 2) + 2*(Points[i+1][1]-Points[i][1])) / (2.0*h)
                    - (M[i+1] - M[i]) * h / 6.0);
        S[i][3] = ((M[i] * pow(Points[i+1][0], 3) - M[i+1]* pow(Points[i][0], 3)) / (6.0 * h)
                + (Points[i+1][0] * Points[i][1] - Points[i][0] * Points[i+1][1]) / h
                - h * (M[i] * Points[i+1][0] - M[i+1] * Points[i][0]) / 6.0);
    }
    return S;
}

void PolyDisp(std::vector<std::vector<double>> S, std::vector<std::vector<double>> points){
    for (int i = 0; i < POINTNUM-1; ++i) {
        std::cout << "S(x)=" << S[i][0] << "x^3+" << S[i][1] << "x^2+" << S[i][2] << "x+" << S[i][3];
        std::cout << ",  x in [" << points[i][0] << ", " << points[i+1][0] << "].\n";
    }
}