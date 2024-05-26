#include "FastFourier.h"

int main(){
    //f1, 2^4
    int n = (int)pow(2, 4);
    std::vector<std::complex<double>> f1_1 = GenerateVec(func1, n);
    std::vector<std::complex<double>> f1_res1 = FFT(f1_1, true);
    std::vector<std::complex<double>> f1_1_IFFT = FFT(f1_res1, false);
    std::cout << "f1结果1： "<< std::endl;
    outputVec(f1_res1);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_1.txt)", f1_1, true);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_res1.txt)", f1_res1, false);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_1_IFFT.txt)", f1_1_IFFT, true);
    //f1, 2^7
    n = (int)pow(2, 7);
    std::vector<std::complex<double>> f1_2 = GenerateVec(func1, n);
    std::vector<std::complex<double>> f1_res2 = FFT(f1_2, true);
    std::vector<std::complex<double>> f1_2_IFFT = FFT(f1_res2, false);
    std::cout << "f1结果2： "<< std::endl;
    outputVec(f1_res2);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_2.txt)", f1_2, true);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_res2.txt)", f1_res2, false);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f1_2_IFFT.txt)", f1_2_IFFT, true);
    //f2, 2^7
    std::vector<std::complex<double>> f2 = GenerateVec(func2, n);
    std::vector<std::complex<double>> f2_res = FFT(f2, true);
    std::vector<std::complex<double>> f2_IFFT = FFT(f2_res, false);
    std::cout << "f2结果： "<< std::endl;
    outputVec(f2_res);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f2.txt)", f2, true);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f2_res.txt)", f2_res, false);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f2_IFFT.txt)", f2_IFFT, true);
    //取f2频率域前25%
    std::vector<std::complex<double>> f2_25pct;
    for(int i = 0; i < n; ++i){
        if(i < 0.25 * n)
            f2_25pct.emplace_back(f2_res[i]);
        else
            f2_25pct.emplace_back((0.0, 0.0));
    }
    std::vector<std::complex<double>> f2_25pct_IFFT = FFT(f2_25pct, false);
    saveRes(R"(D:\Projects\ComputationalMethod\hw6\res\f2_25pct.txt)", f2_25pct_IFFT, true);

    return 0;
}