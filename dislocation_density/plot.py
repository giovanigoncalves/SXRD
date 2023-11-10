import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 


df = pd.read_csv("dislocation_density_data_FCC.txt", sep=",")

fig, ax = plt.subplots()
ax.errorbar(x=df["names"], y=df["rho"], yerr=df["std"], capsize=5)


fig1, ax1 = plt.subplots()
ax1.errorbar(df["names"], df["D"], yerr=df["std_grain"], capsize=5, linestyle="")

plt.show()