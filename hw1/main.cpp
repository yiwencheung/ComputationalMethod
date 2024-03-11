#include <iostream>
#include "mirrorCup.hpp"

int main(){
    //P为观察点，Q为物点
    point P={0,0}, Q, T, R;

    //初始化点坐标
    std::cout << "请输入观察点横坐标：";
    std::cin >> P.x;
    std::cout << "请输入物点坐标：";
    std::cin >> Q.x >> Q.y;

    //寻找切点T与像点R
    T = findT(P,Q);
    R = findR(T,Q);

    //输出结果
    std::cout << "反射点T的坐标为:(" << T.x <<"," <<T.y <<")\n";
    std::cout << "像点R的坐标为:(" << R.x <<"," <<R.y <<")" <<std::endl;
}