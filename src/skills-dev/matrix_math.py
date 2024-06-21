"""Throwaway script to explore syntax for matrix multipliatoin with numpy."""

import numpy as np

print('\n\n')

M = 3

N = M + 1

FLOAT64 = np.float64

# If size is an int, the result is 1D -> 2D requires a tuple.
a = np.zeros(shape=[M, N], dtype=FLOAT64)

# Matrix index is zero based.
a[0, 0] = 1
a[1, 1] = 2
a[2, 2] = 3

print(np.shape(a))

b = np.array([1]*M, dtype=FLOAT64).reshape([3, 1])

print(np.shape(b))

c = np.dot(a, b.transpose())

print(c)
