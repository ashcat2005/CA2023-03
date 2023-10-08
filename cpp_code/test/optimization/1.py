import numpy as np

# create a numpy array of floats
arr = np.array([1.23456789, 2.34567891, 3.45678912])

# set the print options to use scientific notation
np.set_printoptions(precision=2, suppress=True, formatter={'float_kind': '{:0.2e}'.format})

# print the mean of the array
print(np.format_float_scientific(np.mean(arr)))