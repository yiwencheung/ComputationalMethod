#include "FastFourier.h"

double func1(double t){
    double res;
    res = 0.7 * sin(2 * PI * 2 * t) + sin(2 * PI * 5 * t);
    return res;
}

double func2(double t){
    double res;
    std::default_random_engine rand_e;
    std::uniform_real_distribution<double> rand_uni(0, 1.0);

    rand_e.seed(t);
    res = 0.7 * sin(2 * PI * 2 * t) + sin(2 * PI * 5 * t) + rand_uni(rand_e);
    return res;
}

std::vector<std::complex<double>> GenerateVec(FunctionPtr func, int n){
    std::vector<std::complex<double>> res;
    for(int i = 0; i < n; i++)
        res.emplace_back(func(static_cast<double> (i) / n));
    return res;
}

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& f, bool NonInverse){
    int size = (int)f.size();
    if(size <= 1)
        return f;

    std::complex<double> wn;
    std::complex<double> w(1,0);
    std::vector<std::complex<double>> g(size, (0,0));
    std::vector<std::complex<double>> f1,f2;
    for(int i = 0; i < size/2; ++i){
        f1.push_back(f[2 * i]);
        f2.push_back(f[2 * i + 1]);
    }
    wn = NonInverse ? exp(-2*PI / size * std::complex<double>(0.0,1.0))
            : exp(2*PI / size * std::complex<double>(0.0,1.0));
    //分治
    std::vector<std::complex<double>> g1 = FFT(f1, NonInverse);
    std::vector<std::complex<double>> g2 = FFT(f2, NonInverse);
    //结果
    for(int k = 0; k < size/2; ++k){
        if(NonInverse){
            g[k] = (g1[k] + w*g2[k]) / std::complex<double>(2,0);
            g[k + size/2] = (g1[k] - w*g2[k]) / std::complex<double>(2,0);
        }
        else{
            g[k] = (g1[k] + w*g2[k]);
            g[k + size/2] = (g1[k] - w*g2[k]);
        }
        w *= wn;
    }
    return g;
}

void outputVec(const std::vector<std::complex<double>> vec){
    int i = 0;
    for(auto item: vec){
        std::cout << "实部：" << std::real(item) << " 虚部：" << std::imag(item) << "\t";
        i++;
        if(i == 4) {
            std::cout << std::endl;
            i = 0;
        }
    }
}

void saveRes(std::string path, std::vector<std::complex<double>> vec, bool saveReal){
    std::ofstream fout(path, std::ios::out);
    if(!fout.is_open()){
        std::cout << "open file failed!";
        exit(255);
    }
    for(auto item: vec){
        if(saveReal)
            fout << std::real(item) << " " << std::endl;
        else
            fout << std::abs(item) << " " << std::endl;
    }
}