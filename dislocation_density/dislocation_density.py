import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from sklearn.linear_model import LinearRegression
import os



names = ["idx",
         "A111_angle", "A111_area", "A111_fwhm",
         "M101_angle", "M101_area", "M101_fwhm",	
         "M110_angle", "M110_area", "M110_fwhm",	
         "A200_angle", "A200_area", "A200_fwhm",
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

n=200 # Number of the diffraction data
fwhm = [df.A111_fwhm[n], df.A200_fwhm[n], df.A220_fwhm[n], df.A222_fwhm[n]]
fwhm = np.array([float(i) for i in fwhm])

angle = [df.A111_angle[n], df.A200_angle[n], df.A220_angle[n], df.A222_angle[n]]
angle = np.array([np.deg2rad(i/2) for i in angle])


class DislocationDensityMWH:
  
    def __init__(self):
        self.wavelength = 1.4235 # in nanometers
        self.planes = ["111", "200", "220", "222"]
    
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


    def H_squared(self):
        
        h_values = {}
        for plane in self.planes:
            indices = list(plane)
            indices = [int(i) for i in indices]
            h = indices[0]
            k = indices[1]
            l = indices[2]
            numerador = h**2 * k**2 + h**2 * l**2 + k**2 * l**2
            denominador = (h**2 + k**2 + l**2)**2
            h_values[plane] = numerador/denominador
            # print(f"plane: {plane}\n value: {h_values[plane]} ")
        return h_values


    def alpha_value(self, angle, fwhm):
        alpha = np.arange(.0, .01, .000001) # values of alpha to try
        k = self.K_values(angle, fwhm) # it is necessary the delta_k and k values to calculate the y values
        h_squared = self.H_squared().values() # x values for the plot
        x = [i for i in h_squared]
        print(x)
        best_alpha = 0
        best_r_sq = 0
        coef = 0
        intercept = 0
        for num in alpha:
            y = (k["delta_k"]**2 - num) / k["k"]**2
            x = np.array(x).reshape((-1, 1))
            model = LinearRegression().fit(x, y)
            r_sq = model.score(x, y)
            print(r_sq)
            if r_sq > best_r_sq:
                best_r_sq = r_sq
                best_alpha = num
                coef = model.coef_
                intercept = model.intercept_
        print(f"Best value: {best_alpha}")
        print(f"Best r_sq: {best_r_sq}")
        return dict(alpha=best_alpha, r_sq=best_r_sq, coef=coef, intercept=intercept)    
    
    
    def y_values(self, angle, fwhm):
        k = self.K_values(angle, fwhm)
        alpha = self.alpha_value(angle, fwhm)
        return dict(y=(k["delta_k"]**2 - alpha["alpha"]) / k["k"]**2, 
                    alpha=alpha["alpha"], 
                    r_sq=alpha["r_sq"], 
                    coef=alpha["coef"],
                    intercept=alpha["intercept"])
        
        
        
a = DislocationDensityMWH()
k = a.K_values(angle, fwhm)
h_squared = a.H_squared()
h_squared = [value for value in h_squared.values()]
alpha = a.alpha_value(angle, fwhm)
y_values = a.y_values(angle, fwhm)


print(h_squared)


planes = ["111", "200", "220", "222"]

fig, ax = plt.subplots()
ax.set_title("FCC", size=15)
ax.scatter(k["k"], k["delta_k"], facecolors="none", edgecolors="k", marker="*", s=95)
ax.set_xlabel(r"K (nm$^{-1}$)", size=15)
ax.set_ylabel(r"$\Delta$K (nm$^{-1}$)", size=15)
ax.tick_params(axis="both", labelsize=12)
ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)

for i in range(len(planes)):
    ax.annotate("{" + planes[i] + "}", (1.01*k["k"][i], 1.01*k["delta_k"][i]), fontsize=12)
ax.set_xlim(0, .12)
ax.set_ylim(0., .12)

fig.tight_layout()




fig1, ax1 = plt.subplots()
ax1.set_title("FCC", size=15)
ax1.scatter(h_squared, y_values["y"], facecolors="none", edgecolors="k", marker="*", s=95)
ax1.set_xlabel(r"H$^2$", size=15)
ax1.set_ylabel(r"$\frac{\Delta K^2 - \alpha}{k^2}$", size=15)
ax1.tick_params(axis="both", labelsize=12)
ax1.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)

for i in range(len(planes)):
    ax1.annotate("{" + planes[i] + "}", (1.01*h_squared[i], 1.01*y_values["y"][i]), fontsize=12)
    
ax1.annotate(f"R²: {y_values['r_sq']:.4f}", (.25, y_values["y"].max()), fontsize=12)
ax1.annotate(rf"$\alpha$: {y_values['alpha']}", (.25, y_values["y"].max()*.975), fontsize=12)

# ax.set_xlim(0, .12)
# ax.set_ylim(0., .2)

x1 = np.arange(-.02, .35, .01)
ax1.plot(x1, y_values["coef"] * x1 + y_values["intercept"])
ax1.annotate(rf"q: {- float(y_values['coef']) / float(y_values['intercept']):.3f}", (.25, y_values["y"].max()*.95), fontsize=12)


fig1.tight_layout()

plt.show()