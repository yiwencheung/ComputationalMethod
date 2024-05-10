#include "Spline.h"

int main(){
    //计算三次样条插值
    int h = 1;      //步长
    int i = 0;
    std::fstream fin;
    std::string rdline;
    std::vector<std::vector<double>> points(POINTNUM, std::vector<double>(DIV, 0));
    //读取并分割数据
    fin.open("point.txt", std::ios::in);
    while(std::getline(fin, rdline)){
        std::vector<std::string> data = split(rdline, " ");
        points[i][0] = std::stod(data[0]);
        points[i][1] = std::stod(data[1]);
        ++i;
    }
    


    return 0;
}