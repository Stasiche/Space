import numpy as np
import matplotlib.pyplot as plt
import random

TOTAL = 5
STEP = 2


def func(x):
    return x - 1


def func2(x):
    return 2*x + 1


def func3(x):
    return 1.2*x + 0.2


X1 = [i for i in range(0, 5, 2)]
Y1 = func(np.array(X1))

X2 = [i for i in range(1, 5, 2)]
Y2 = func2(np.array(X2))

X_real = np.arange(0, 5, 0.001)
Y_real = func(np.array(X1))
ax = plt.axis([-2, 5, -3, 8])

A = np.empty((TOTAL, 2))
A[:, 0] = 1
A[:, 1] = np.array(X1+X2)

theta = np.linalg.pinv(A).dot(np.array(list(Y1)+list(Y2)))
Y_prediction = A.dot(theta)

print(theta)

plt.plot(X1+X2, list(Y1)+list(Y2), 'bo')
plt.grid()
plt.plot(np.arange(-10, 5, 0.001), func(np.arange(-10, 5, 0.001)), 'g', linewidth=2.0)
plt.plot(np.arange(-10, 5, 0.001), func2(np.arange(-10, 5, 0.001)), 'r', linewidth=2.0)
plt.plot(np.arange(-10, 5, 0.001), func3(np.arange(-10, 5, 0.001)), 'b', linewidth=2.0)

plt.plot(X1+X2, func3(np.array(X1+X2)), 'bo')

plt.show()
