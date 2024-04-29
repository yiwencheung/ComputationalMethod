import numpy as np
import matplotlib.pyplot as plt

FILE_PATH = "D:\\Projects\\ComputationalMethod\\hw4\\res\\"

def main():
    # 读取文件数据
    with open(FILE_PATH + "label0.txt", 'r') as f0:
        data0 = [line.strip().split(',') for line in f0]
        x0 = [float(row[0]) for row in data0]
        y0 = [float(row[1]) for row in data0]
    with open(FILE_PATH + "label1.txt", 'r') as f1:
        data1 = [line.strip().split(',') for line in f1]
        x1 = [float(row[0]) for row in data1]
        y1 = [float(row[1]) for row in data1]
    with open(FILE_PATH + "label2.txt", 'r') as f2:
        data2 = [line.strip().split(',') for line in f2]
        x2 = [float(row[0]) for row in data2]
        y2 = [float(row[1]) for row in data2]

    #可视化处理
    plt.scatter(x0, y0, color = "red" ,label = "label = 0")
    plt.scatter(x1, y1, color = "blue", label = "label = 1")
    plt.scatter(x2, y2, color = "green", label = "label = 2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Projection of datas in the principal direction")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()



