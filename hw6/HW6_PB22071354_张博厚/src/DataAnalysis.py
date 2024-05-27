import numpy as np
import matplotlib.pyplot as plt

ROOT_PATH = "D:\\Projects\\ComputationalMethod\\hw6\\res\\datas\\"
SAVE_PATH = "D:\\Projects\\ComputationalMethod\\hw6\\res\\figs\\"


def txt2array(path):
    with open(path, 'r') as f:
        data = [float(num) for num in f.read().split()]
        data_array = np.array(data)
    return data_array


def plot_FFT(title, y, fre):
    x = np.array([i for i in range(0, fre)])
    plt.plot(x, y, label=f"frequency = {fre}", linewidth=2)
    plt.xlabel("fre")
    plt.ylabel("|g_i|")
    plt.title(title)
    plt.legend()
    plt.savefig(SAVE_PATH + title + ".png")
    plt.close()


def plot_IFFT(title, y1, y2, fre):
    x = np.array([i / fre for i in range(0, fre)])
    plt.plot(x, y1, label='origin', color='red', linewidth=2)
    plt.plot(x, y2, label='IFFT', color='blue', linewidth=2)
    plt.xlabel("t")
    plt.ylabel("real")
    plt.title(title)
    plt.legend()
    plt.savefig(SAVE_PATH + title + ".png")
    plt.close()


def plot_IFFT_25(title, y1, y2, y3, fre):
    x = np.array([i / fre for i in range(0, fre)])
    plt.plot(x, y1, label='origin', color='red', linewidth=2)
    plt.plot(x, y2, label='IFFT', color='blue', linewidth=2)
    plt.plot(x, y3, label='IFFT_0.25', color='green', linewidth=2)
    plt.xlabel("t")
    plt.ylabel("real")
    plt.title(title)
    plt.legend()
    plt.savefig(SAVE_PATH + title + ".png")
    plt.close()


def main():
    #读取数据
    f1_origin1 = txt2array(ROOT_PATH + "f1_1.txt")
    f1_origin2 = txt2array(ROOT_PATH + "f1_2.txt")
    f1_FFT1 = txt2array(ROOT_PATH + "f1_res1.txt")
    f1_FFT2 = txt2array(ROOT_PATH + "f1_res2.txt")
    f1_IFFT1 =  txt2array(ROOT_PATH + "f1_1_IFFT.txt")
    f1_IFFT2 = txt2array(ROOT_PATH + "f1_2_IFFT.txt")
    f2_origin = txt2array(ROOT_PATH + "f2.txt")
    f2_FFT = txt2array(ROOT_PATH + "f2_res.txt")
    f2_IFFT = txt2array(ROOT_PATH + "f2_IFFT.txt")
    f2_25pct = txt2array(ROOT_PATH + "f2_25pct.txt")
    #画图
    plot_FFT("FFT_f1_n1", f1_FFT1, 16)
    plot_FFT("FFT_f1_n2", f1_FFT2, 128)
    plot_FFT("FFT_f2", f2_FFT, 128)
    plot_IFFT("IFFT_f1_n1", f1_origin1, f1_IFFT1, 16)
    plot_IFFT("IFFT_f1_n2", f1_origin2, f1_IFFT2, 128)
    plot_IFFT_25("IFFT_f2", f2_origin, f2_IFFT, f2_25pct, 128)


if __name__ == '__main__':
    main()