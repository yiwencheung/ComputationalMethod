import numpy as np
import matplotlib.pyplot as plt

FILE_PATH = "D:\\Projects\\ComputationalMethod\\hw5\\res\\"

# 从文件中读取参数
def read_parameters(filename):
    with open(filename, 'r') as file:
        parameters = [[float(num) for num in line.split()] for line in file]
    return parameters

# 定义20个三次函数
def cubic_function(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d

def piecewise_cubic(x, parameters):
    conditions = []
    functions = []
    pos = -9
    for a, b, c, d in parameters:
        conditions.append((x >= pos) & (x < pos+1))
        functions.append(lambda x, a=a, b=b, c=c, d=d: cubic_function(x, a, b, c, d))
        pos += 1
    return np.piecewise(x, conditions, functions)


def main():
    parameters_0 = read_parameters(FILE_PATH + "para0.txt")
    parameters_1 = read_parameters(FILE_PATH + "para1.txt")
    x = np.linspace(-9,11,1000)
    y0 = piecewise_cubic(x, parameters_0)
    y1 = piecewise_cubic(x, parameters_1)
    # 绘制图形
    plt.plot(x, y0, 'b', label = 'before')
    plt.plot(x, y1, 'r', label = 'after')
    plt.title('Comparison')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
