import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import os

planes = ["110", "200", "211", "220"]

names = ["idx",
         "A111_angle_", "A111_area", "A111_fwhm",
         "M101_angle", "M101_area", "M101_fwhm",	
         "M110_angle", "M110_area", "M110_fwhm",	
         "A002_angle", "A002_area", "A002_fwhm",
         "M002_angle", "M002_area", "M002_fwhm",	
         "M200_angle", "M200_area", "M200_fwhm",	
         "A220_angle", "A220_area", "A220_fwhm",	
         "M112_angle", "M112_area", "M112_fwhm",	
         "M211_angle", "M211_area", "M211_fwhm",	
         "A311_angle", "A311_area", "A311_fwhm",	
         "A222_angle", "A222_area", "A222_fwhm",	
         "M202_angle", "M202_area", "M202_fwhm",	
         "M220_angle", "M220_area", "M220_fwhm",	
        ]

df = pd.read_excel("gg_qp_250_2h.xlsx", skiprows=1, names=names)
df.replace("NP", np.nan, inplace=True)

fwhm = [df.M110_fwhm[0], df.M200_fwhm[0], df.M211_fwhm[0], df.M220_fwhm[0]]
fwhm = np.array([float(i) for i in fwhm])

angle = [df.M110_angle[0], df.M200_angle[0], df.M211_angle[0], df.M220_angle[0]]
angle = np.array([np.deg2rad(i/2) for i in angle])


class DislocationDensityMWH:
  
    def __init__(self):
        self.wavelength = 1.4235 # in nanometers
    
    
    def dislocation_density(D, A, b, p, K, C, A_prime, Q, beta_prime, Wg):
        
        x = 0.9/D
        y = np.sqrt(np.pi * A**2 * b**2 / 2) * np.sqrt(p) * (K * np.sqrt(C))
        z = (np.pi * A_prime * b**2 / 2) * np.sqrt(Q) * (K**2 * C)
        w = beta_prime * Wg
        
        return dict("Dislocation Density", x + y + z + w)


    def K_values(self, angle, fwhm):
        """
            Calculates the K and Delta_K values. 
            Delta_K values are reffered to FWHM data.
            Return a dict with both values.
        """       
        k = 2 * np.sin(angle) / (self.wavelength)
        delta_k = 2 * np.cos(angle) * fwhm / self.wavelength
        
        return dict(k=k, delta_k=delta_k)

    def H_squared(self, plane):
        indices = list(plane)
        indices = [int(i) for i in indices]
        h = indices[0]
        k = indices[1]
        l = indices[2]
        numerador = h**2 * k**2 + h**2 * l**2 + k**2 * l**2
        denominador = (h**2 + k**2 + l**2)**2
        
        return numerador/denominador

a = DislocationDensityMWH()
k = a.K_values(angle, fwhm)


fig, ax = plt.subplots()
ax.scatter(k["k"], k["delta_k"], facecolors="none", edgecolors="k", marker="*", s=95)
ax.set_xlabel(r"K (nm$^{-1}$)", size=15)
ax.set_ylabel(r"$\Delta$K (nm$^{-1}$)", size=15)
ax.tick_params(axis="both", labelsize=12)

for i in range(len(planes)):
    ax.annotate("{" + planes[i] + "}", (1.01*k["k"][i], 1.01*k["delta_k"][i]), fontsize=12)

fig.tight_layout()

plt.show()