import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 


df = pd.read_csv("dislocation_density_data_BCC.txt", sep=",")

fig, ax = plt.subplots()
ax.errorbar(x=df["names"], y=df["rho"], yerr=df["std"], 
            capsize=5, 
            ls="--", 
            c="k")
ax.set_xlabel("Measurement", size=15)
ax.set_ylabel(r"Dislocation Density (m$^{-2}$)", size=15)
ax.set_yscale("log")
ax.tick_params(axis="both", labelsize=12)
fig.tight_layout()

fig1, ax1 = plt.subplots()
ax1.errorbar(df["names"], df["D"]/1000, yerr=df["std_grain"]/1000, 
             capsize=5, 
             linestyle="--", 
             c="k")
ax1.set_xlabel("Measurement", size=15)
ax1.tick_params(axis="both", labelsize=12)
ax1.set_ylabel(r"Crystallite size ($\mu$m)", size=15)
fig1.tight_layout()

plt.show()