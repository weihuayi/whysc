import numpy as np
import matplotlib.pyplot as plt

def Histogram_plot(b):
    a = np.abs(np.array(b)-1)
    A = np.linspace(0, 1, 51)
    B = np.zeros(50)
    for i in range(50):
        flag = (a<= A[i+1]) & (a>A[i])
        B[i] = flag.sum()
        s = 0
    plt.bar(A[:-1], B, 0.018)
    plt.ylim(0, 70)
    plt.show()



