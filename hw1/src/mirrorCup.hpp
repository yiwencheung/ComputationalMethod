#include <iostream>
#include <cmath>

typedef struct{
    double x;
    double y;
} point;

double computeDist(point A, point B);
point findT(point P, point Q);
point findR(point T, point Q);

