#include <iostream>
#include "mirrorCup.hpp"

int main(){
    //PΪ�۲�㣬QΪ���
    point P={0,0}, Q, T, R;

    //��ʼ��������
    std::cout << "������۲������꣺";
    std::cin >> P.x;
    std::cout << "������������꣺";
    std::cin >> Q.x >> Q.y;

    //Ѱ���е�T�����R
    T = findT(P,Q);
    R = findR(T,Q);

    //������
    std::cout << "�����T������Ϊ:(" << T.x <<"," <<T.y <<")\n";
    std::cout << "���R������Ϊ:(" << R.x <<"," <<R.y <<")" <<std::endl;
}