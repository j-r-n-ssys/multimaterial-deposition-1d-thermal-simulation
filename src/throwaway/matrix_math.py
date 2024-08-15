"""Throwaway script to explore syntax for matrix multipliatoin with numpy."""

import numpy as np

# num rows
M = 3

# num col
N = M + 1

FLOAT64 = np.float64

# If size is an int, the result is 1D -> 2D requires a tuple.
# Initialize an M x N matrix (conduction matrix). 
a = np.zeros(shape=[M, N], dtype=FLOAT64)

# Matrix index is zero based.
a[0, 0] = 1
a[1, 1] = 2
a[2, 2] = 3

# Initialize an N x 1 vector (future: temperature in Z).
b = np.array([1]*N, dtype=FLOAT64)

print(f'shape of a is {a.shape}')

print(f'shape of b is {b.shape}')

c = np.dot(a, b.T)

print(c)


m = np.array([1,2,3], dtype=FLOAT64)

print((m * m))


print('\n\n')

a = np.array([1,2,3], dtype=FLOAT64)

b = np.array([1,2,3], dtype=FLOAT64)

print(a)
print(b)
print(np.divide(a,b))


d = np.array([1,2,3])



e = np.ones(shape = d.shape)

print(np.shape(d))

print(np.shape(e))

e = np.divide(e,d)


print(e)

f = np.divide(1, d)

print(f)

g = np.array([1,2,3])

print(np.power(10, g))