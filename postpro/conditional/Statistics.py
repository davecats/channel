# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:23:31 2022

@author: Test
"""

import numpy as np
import matplotlib.pyplot as plt


def read_file(filename):

    dtype_a = np.dtype(np.intc)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        # read the date into an array of double
        data = np.fromfile(binary_file, dtype_a)

        return data


# writes the zero-crossing positions for both positive to negative and vice-versa
array_pos = np.zeros([], dtype=np.intc)

array_neg = np.zeros([], dtype=np.intc)

array_pos = read_file("pos_neg.out")

array_neg = read_file("neg_pos.out")

array = np.zeros([], dtype=np.intc)

array = np.append(array_pos, array_neg)

array = np.sort(array, axis=None)



print('Kindly find the combined zero-crossing positions below: \n')

print(array)
print('\n')

diff_array = np.zeros([len(array)-1], dtype=np.intc)

i = 1
while i < len(array):
    diff_array[i-1] = array[i] - array[i-1]
    i += 1
    
print('Kindly find the interval between succesive zero-crossing positions below: \n')

print(diff_array)
print('\n') 

max_diff = np.max(diff_array)
print('The maximum interval is ' + str(max_diff) + '\n')

min_diff = np.min(diff_array)
print('The minimum interval is ' + str(min_diff) + '\n')

average_diff = np.mean(diff_array) 
print('The average of the intervals is ' + str(round(average_diff, 3)) + '\n')

median_diff = np.median(diff_array)
print('The median of the intervals is ' + str(median_diff) + '\n')

First_quartile_diff = np.quantile(diff_array, 0.25)
print('The 1st quartile of the intervals is ' + str(First_quartile_diff) + '\n')

Third_quartile_diff = np.quantile(diff_array, 0.75)
print('The 3rd quartile of the intervals is ' + str(Third_quartile_diff) + '\n')

Int_Range_diff = Third_quartile_diff - First_quartile_diff

print('The interquartile range of the intervals is ' + str(Int_Range_diff) + '\n')













    