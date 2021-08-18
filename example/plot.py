import numpy as np
import matplotlib.pyplot as plt

def Histogram_plot(b):
    a = np.array(b)
    a = 1/a
    L = len(b)
    A = np.linspace(0, 1, 101)
    B = np.zeros(100)
    for i in range(100):
        flag = (a<= A[i+1]) & (a>A[i])
        B[i] = flag.sum()/L
        s = 0
    plt.bar(A[:-1], B, 0.009)
    #plt.ylim(0, np.max(B)*1.05)
    plt.ylim(0, 0.1)
    plt.xlabel("cell quality")
    plt.ylabel("Proportion")
    plt.show()

