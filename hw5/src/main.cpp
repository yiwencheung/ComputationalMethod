#include "Spline.h"

int main(){
    /*��������������ֵ
     * ʹ�ô�M����ȡ��Ȼ�߽�������
    */
    //���õ��ı���
    int h = 1;      //����
    int i = 0;
    double lambda = 0.5, mu = 0.5;
    std::fstream fin;
    std::string rdline;
    std::vector<std::vector<double>> points(POINTNUM, std::vector<double>(DIV, 0));
    std::vector<double> M(POINTNUM, 0);
    std::vector<std::vector<double>> S(POINTNUM - 1, std::vector<double>(POWER + 1, 0));
    //��ȡ���ָ�����
    fin.open(FILE_PATH, std::ios::in);
    while(std::getline(fin, rdline)){
        std::vector<std::string> data = split(rdline, "\t");
        points[i][0] = std::stod(data[0]);
        points[i][1] = std::stod(data[1]);
        ++i;
    }
    fin.close();
    //���M
    M = CubicSpline_natural(lambda, mu, h, points);
    //��M������ʽϵ��
    S = SolveS(M, points, h);
    //�������ʽ
    PolyDisp(S, points);

    //�������λ�����¼���
    std::vector<double> M_plus(POINTNUM, 0);
    std::vector<std::vector<double>> S_plus(POINTNUM - 1, std::vector<double>(POWER + 1, 0));
    points[9][1] = 10;
    M_plus = CubicSpline_natural(lambda, mu, h, points);
    S_plus = SolveS(M_plus, points, h);
    std::cout << "\n���ĵ�ʮ����󣬽��Ϊ:" << std::endl;
    PolyDisp(S_plus, points);

    //�Ƚ�ǰ�����ν��
    std::cout << "����ǰ�󣬶���ʽ��ϵ���仯Ϊ��\n";
    for(int j = 0; j < S.size(); ++j){
        std::cout << S[j][0] - S_plus[j][0] << ", ";
        std::cout << S[j][1] - S_plus[j][1] << ", ";
        std::cout << S[j][2] - S_plus[j][2] << ", ";
        std::cout << S[j][3] - S_plus[j][3] << std::endl;
    }

    //������
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