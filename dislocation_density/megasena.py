import numpy as np
import pandas as pd 
from random import choice
from pprint import pprint

options = np.arange(1,61)

# print(options)


sorted_numbers = {}
for i in range(6):
    option = list(options)
    l = []
    for j in range(6):
        number = choice(option)
        l.append(number)
        option.remove(number)
    l.sort()
    sorted_numbers[i] = l      
        
pprint(sorted_numbers)