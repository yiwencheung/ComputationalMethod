#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>

#define PI 3.1415926535897932

typedef double (*FunctionPtr)(double);

//本实验所用到的两个函数
double func1(double t);
double func2(double t);
//输入函数与n，得到f数组
std::vector<std::complex<double>> GenerateVec(FunctionPtr func, int n);
/*快速傅里叶变换
 * NonInverse为真表示FFT， 为假表示IFFT
 */
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& f, bool NonInverse);

void outputVec(const std::vector<std::complex<double>> vec);
void saveRes(std::string path, std::vector<std::complex<double>> vec, bool saveReal);