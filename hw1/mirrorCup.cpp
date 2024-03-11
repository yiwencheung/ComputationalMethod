#include <iostream>
#include "mirrorCup.hpp"

point findT(point P, point Q){
    point T={0,0};
    //首先寻找角POT的上限
    double POG , POQ;
    POG = acos(1/-P.x);
    POQ = acos(-Q.x/ sqrt(Q.x*Q.x + Q.y*Q.y));
    //用二分法寻找点T
    double POT;
    double low=0, high = fmin(POG,POQ);
    while(high-low > 1e-15){
        //取中间角度并计算T
        POT = (low+high)/2;
        T.x = -cos(POT);
        T.y = sin(POT);
        //计算角PTV，QTV并比较
        double l_PT = computeDist(P,T), l_QT = computeDist(Q,T);
        double PTV = asin(fabs(T.y*P.x-T.x*P.y) / l_PT);
        double QTV = asin(fabs(T.y*Q.x-T.x*Q.y) / l_QT);
        if(fabs(PTV-QTV)<1e-15)
            break;
        else if(PTV > QTV)
            high = POT;
        else
            low = POT;
    }
    return T;
}

point findR(point T, point Q){
    //R为像点， M为对称点
    point R, M;
    M.x = T.x * T.y * (T.y-Q.y) + T.x * T.x * T.x + T.y * T.y * Q.x;
    M.y = T.y - T.x * (M.x - T.x) /T.y;
    //利用对称性求解R
    R.x = 2*M.x - Q.x;
    R.y = 2*M.y - Q.y;

    return R;
}

double computeDist(point A, point B){
    return sqrt((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y));
}