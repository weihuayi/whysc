import numpy as np
import time
import matplotlib.pyplot as plt

def Histogram_plot(b):
    print("最差单元质量：", np.max(b))
    print("最优单元质量：", np.min(b))
    print("平均单元质量：", np.average(b))
    a = np.array(b)
    a = 1/a
    L = len(b)
    A = np.linspace(0, 1, 51)
    B = np.zeros(50)
    for i in range(50):
        flag = (a<= A[i+1]) & (a>A[i])
        B[i] = flag.sum()/L
        s = 0
    fig = plt.figure()
    axes = fig.gca()
    axes.bar(A[:-1], B, 0.018)
    #axes.ylim(0, np.max(B)*1.05)
    plt.ylim(0, 0.1)
    plt.xlabel("cell quality")
    plt.ylabel("Proportion")
    fig.savefig('%f.svg' %B[25], dpi=600)
    plt.show()

