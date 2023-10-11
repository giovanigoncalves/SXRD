import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from sklearn.linear_model import LinearRegression
import os
from pprint import pprint

# Name of coluns for data containing the FWHM values
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

# Loading the data containing FWHM values
df = pd.read_excel("gg_qp_250_2h.xlsx", skiprows=1, names=names)
df.replace("NP", np.nan, inplace=True)


# Chosing one unique diffraction spectra to calc dislocation density (FWHM and angle info)
SPECTRAS = np.arange(900, 1000, 100) # Number of the diffraction data


fwhm = {}
angle = {}
angle_1 = {}
angle_2 = {}
for i in SPECTRAS:
    fwhm[i] = [df.M110_fwhm[i], df.M200_fwhm[i], df.M211_fwhm[i], df.M220_fwhm[i]]
    fwhm[i] = np.array([float(n) for n in fwhm[i]])

    angle[i] = [df.M110_angle[i], df.M200_angle[i], df.M211_angle[i], df.M220_angle[i]]
    angle[i] = np.array([np.deg2rad(n/2) for n in angle[i]]) # Transforming the angles from degree to radians

    angle_1[i] = [df.M110_angle[i] - df.M110_fwhm[i] / 2, 
                  df.M200_angle[i] - df.M200_fwhm[i] / 2, 
                  df.M211_angle[i] - df.M211_fwhm[i] / 2, 
                  df.M220_angle[i] - df.M220_fwhm[i] / 2
                  ]
    angle_1[i] = np.array([np.deg2rad(n/2) for n in angle_1[i]]) # Transforming the angles from degree to radians
    angle_2[i] = [df.M110_angle[i] + df.M110_fwhm[i] / 2, 
                  df.M200_angle[i] + df.M200_fwhm[i] / 2, 
                  df.M211_angle[i] + df.M211_fwhm[i] / 2, 
                  df.M220_angle[i] + df.M220_fwhm[i] / 2
                  ]
    angle_2[i] = np.array([np.deg2rad(n/2) for n in angle_2[i]]) # Transforming the angles from degree to radians
    

# Here we start the play

class DislocationDensityMWH:
    # Initial data (wavelength and phase planes for BCC)
    def __init__(self):
        self.wavelength = 1.4235 # in nanometers
        self.planes = ["110", "200", "211", "220"]
        self.a_bcc = 2.87
        self.a_fcc = 3.59
        self.alpha_range = np.arange(0, 1e-3, 1e-8)
    # Calculating the K values for the different planes
    def k_values(self, angle, angle_1, angle_2):
        """
            Calculates the K and Delta_K values. 
        """
        k = {}
        delta_k = {}
        for key, value in angle.items():     
            k[key] = 2 * np.sin(value) / (self.wavelength) # the angles need to be in radians
            
        for i in SPECTRAS:
            delta_k[i] = 2 * (np.sin(angle_2[i] - np.sin(angle_1[i]) / self.wavelength))
            

        
        k_values = pd.DataFrame()
        k_values["measurement"] = k.keys()

        for measure in range(len(self.planes)):
            values = []
            for i in k:
                values.append(k[i][measure])
            k_values[self.planes[measure]] = values
        k_values.set_index("measurement", inplace=True)

        delta_k_values = pd.DataFrame()
        delta_k_values["measurement"] = delta_k.keys()
        for measure in range(len(self.planes)):
            values = []
            for i in delta_k:
                values.append(delta_k[i][measure])
            delta_k_values[self.planes[measure]] = values
        delta_k_values.set_index("measurement", inplace=True)

        return k_values, delta_k_values
    

    def H_squared(self):
        
        h_values = {}
        for plane in self.planes:
            indices = list(plane)
            indices = [int(i) for i in indices]
            h = indices[0]
            k = indices[1]
            l = indices[2]
            h_values[plane] = ( h**2 * k**2 + h**2 * l**2 + k**2 * l**2 ) / ( h**2 + k**2 + l**2 ) **2
            # print(f"plane: {plane}\n value: {h_values[plane]} ")
        return h_values


    def alpha_value(self, angle, angle_1, angle_2):
        k, delta_k = self.k_values(angle, angle_1, angle_2)
        
        alpha = self.alpha_range # values of alpha to try
        h_squared = self.H_squared().values() # x values for the plot
        new_x = np.arange(0, 1, .01)
        x = [i for i in h_squared]
        x = np.array(x).reshape((-1, 1))

        # print(x)
        best_alpha = {}
        best_r_sq = {}
        coef = {}
        intercept = {}
        q = {}
        
        delta_k = delta_k.T
        k = k.T
        
        for col in delta_k.columns:
            p = 0
            best_alpha[col] = 0
            best_r_sq[col] = 0
            coef[col] = 0
            intercept[col] = 0
            for num in alpha:
                y =  (delta_k[col]**2 - num) / k[col]**2
                # print(y)
                model = LinearRegression().fit(x, y)
                r_sq = model.score(x, y)
                # print(r_sq)
                if r_sq > best_r_sq[col]:
                    best_r_sq[col] = r_sq
                    best_alpha[col] = num
                    coef[col] = model.coef_
                    intercept[col] = model.intercept_
                p += 1
                print(f"{p*100/(len(alpha)*len(delta_k.columns)):.2f} %")
            q[col] = ( -1 * coef[col] ) / intercept[col]            
            
                
            fig, ax = plt.subplots()
            ax.scatter(h_squared,  ( (delta_k[col] - best_alpha[col]) / k[col] )**2, c="k")
            ax.plot(new_x, intercept[col] + coef[col]*new_x, c="k")
            ax.set_xlabel(r"H$^2$", size=15)
            ax.set_ylabel(r"$\frac{(\Delta K^2 - \alpha)}{K^2}$", size=15)
            ax.tick_params(axis="both", labelsize=12)
            ax.set_ylim(0, )
            ax.annotate(f"Measure: {col}", (0, .01))
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
            fig.tight_layout()
            
            try:
                os.mkdir("q_results_plot")
            except FileExistsError:
                pass
            
            plt.savefig(f"q_results_plot/Measure_{col}.png")
            
            plt.show()
            
            print(f"Best value for {col}: {best_alpha[col]}")
            print(f"Best r_sq for {col}: {best_r_sq[col]:.2f}")
            print(f"1/q for {col}: {-1 * intercept[col] / coef[col]}")
            print(f"q for {col}: {-1 * coef[col] / intercept[col]}")
            
        return dict(alpha=best_alpha, r_sq=best_r_sq, coef=coef, intercept=intercept, q=q)    
    
    
    def q_values(self, angle, angle_1, angle_2, structure="BCC"):
        
        q = self.alpha_value(angle, angle_1, angle_2)["q"]
        
        if structure.upper() == "BCC":
            df_screw = pd.read_csv(f"abcd_parameters_for_q_screw_dislocation_bcc.txt", delimiter=";")
            df_edge = pd.read_csv(f"abcd_parameters_for_q_edge_dislocation_bcc.txt", delimiter=";")
            # c11, c12 and c44 were chosen from F. HajyAkbary et al. (2015) Materials Science & Engineering A 639
            c11 = 230
            c12 = 135
            c44 = 117
            
            if abs(c12/c44 - 0.5) < abs(c12/c44 - 1) < abs(c12/c44 - 2):
                a_q_screw = df_screw["c12/c44=0.5"][0]
                b_q_screw = df_screw["c12/c44=0.5"][1]
                c_q_screw = df_screw["c12/c44=0.5"][2]
                d_q_screw = df_screw["c12/c44=0.5"][3]  
                
                a_q_edge = df_edge["c12/c44=0.5"][0]
                b_q_edge = df_edge["c12/c44=0.5"][1]
                c_q_edge = df_edge["c12/c44=0.5"][2]
                d_q_edge = df_edge["c12/c44=0.5"][3]

            elif abs(c12/c44 - 1) < abs(c12/c44 - 2):
                a_q_screw = df_screw["c12/c44=1"][0]
                b_q_screw = df_screw["c12/c44=1"][1]
                c_q_screw = df_screw["c12/c44=1"][2]
                d_q_screw = df_screw["c12/c44=1"][3]
                
                a_q_edge = df_edge["c12/c44=1"][0]
                b_q_edge = df_edge["c12/c44=1"][1]
                c_q_edge = df_edge["c12/c44=1"][2]
                d_q_edge = df_edge["c12/c44=1"][3]

            else:
                a_q_screw = df_screw["c12/c44=2"][0]
                b_q_screw = df_screw["c12/c44=2"][1]
                c_q_screw = df_screw["c12/c44=2"][2]
                d_q_screw = df_screw["c12/c44=2"][3]
                
                a_q_edge = df_edge["c12/c44=2"][0]
                b_q_edge = df_edge["c12/c44=2"][1]
                c_q_edge = df_edge["c12/c44=2"][2]
                d_q_edge = df_edge["c12/c44=2"][3]

   
        elif structure.upper() == "FCC":
            name = "fcc"
            
        else:
            print("Invalid structure! Chose BCC or FCC!")
                   

        A = 2 * c44 / (c11 - c12) # Elastic anisotropy
        q_screw = a_q_screw * (1 - np.exp(-A/b_q_screw)) + c_q_screw * A + d_q_screw # constant value for 100% screw dislocations
        q_edge = a_q_edge * (1 - np.exp(-A/b_q_edge)) + c_q_edge * A + d_q_edge # constant value for 100% edge dislocations

        f_edge = {}
        f_screw = {}
        for key in q:
            f_edge[key] = ( q_screw - q[key] ) / ( q_screw - q_edge )  
            f_screw[key] = 1- f_edge[key]
        
        return dict(f_edge=f_edge, f_screw=f_screw, q_screw=q_screw, q_edge=q_edge, q=q, A=A)

    
    def ch00(self, structure="BCC"):
        q_theoretical = self.q_values()
        
        if structure.upper() == "BCC":
            df_ch00_edge = pd.read_csv("abcd_parameters_for_ch00_edge_dislocation_bcc.txt", delimiter=";")
            c11 = 230
            c12 = 135
            c44 = 117
            
            # Os valores praticamente independem do c12/c44 para screw dislocations in BCC structure
            a_ch00_screw = 0.174
            b_ch00_screw = 1.9522
            c_ch00_screw = 0.0293
            d_ch00_screw = 0.0662

            if abs(c12/c44 - 0.5) < abs(c12/c44 - 1) < abs(c12/c44 - 2):
                a_ch00_edge = df_ch00_edge["c12/c44=0.5"][0]
                b_ch00_edge = df_ch00_edge["c12/c44=0.5"][1]
                c_ch00_edge = df_ch00_edge["c12/c44=0.5"][2]
                d_ch00_edge = df_ch00_edge["c12/c44=0.5"][3]

            elif abs(c12/c44 - 1) < abs(c12/c44 - 2):
                a_ch00_edge = df_ch00_edge["c12/c44=1"][0]
                b_ch00_edge = df_ch00_edge["c12/c44=1"][1]
                c_ch00_edge = df_ch00_edge["c12/c44=1"][2]
                d_ch00_edge = df_ch00_edge["c12/c44=1"][3]

            else:
                a_ch00_edge = df_ch00_edge["c12/c44=2"][0]
                b_ch00_edge = df_ch00_edge["c12/c44=2"][1]
                c_ch00_edge = df_ch00_edge["c12/c44=2"][2]
                d_ch00_edge = df_ch00_edge["c12/c44=2"][3]

        ch00_edge = a_ch00_edge * ( 1 - np.exp( - self.A / b_ch00_edge) ) + c_ch00_edge * self.A + d_ch00_edge
        ch00_screw = a_ch00_screw * ( 1 - np.exp( - self.A / b_ch00_screw) ) + c_ch00_screw * self.A + d_ch00_screw
        
        ch00 = ( q_theoretical["f_edge"] * ch00_edge ) + ( q_theoretical["f_screw"] * ch00_screw )
        
        return ch00
        
        
    def c(self, q):
        q_theoretical = q_theoretical(q, structure="BCC")
        ch00 = self.ch00(structure="BCC")
        h_values = self.H_squared()
        
        df_c = pd.DataFrame()
        df_c["planes"] = self.planes
        df_c.set_index("planes", inplace=True)
        for i in range(len(ch00)):
            c = []
            for j in self.planes:
                c.append = ch00[i] * ( 1 - [i] * h_values[j])
            df_c[i] = c
        
        return df_c
    
    def b(self, structure="BCC"):

        if structure.upper() == "BCC":
            a = self.a_bcc
        elif structure.upper() == "FCC":
            a = self.a_fcc

        b = {}
        for i in self.planes:
            plane = [int(i) for i in list(i)]           
            b[i] = (a / 2) * np.sqrt(plane[0]**2 + plane[1]**2 + plane[2]**2)

        if structure.upper() == "BCC":
            return b["110"] # According to ref. T. Ungár et al. (1999) pag. 993
        elif structure.upper() == "FCC":
            return b["111"] # According to ref. T. Ungár et al. (1999) pag. 993
        

        
a = DislocationDensityMWH()
# k, delta_k = a.k_values(angle, angle_1, angle_2)
# h_squared = a.H_squared()
# h_squared = [value for value in h_squared.values()]
# alpha = a.alpha_value(angle, angle_1, angle_2)
# q = alpha["q"]
q_values = a.q_values(angle, angle_1, angle_2)
# y_values = a.y_values(angle, fwhm)
pprint(q_values)
print()
print()
# print(k)



# print(new)

# print(alpha)
# print()
# print(delta_k.T)


# a = {"angle": [1, 2, 3, 4, 5, 6, 7, 8, 9], "intensity": [10, 20, 30], "failure": [80, 45, 25, 35]}

# for k, v in a.items():
#     print(v)

# # Plotting condition x delta K
# fig, ax = plt.subplots()

# for col in delta_k.columns:
#     ax.scatter(delta_k.index, delta_k[col], label=f"Medida: {col}", marker="^")
# plt.show()



# planes = ["111", "200", "220", "222"]

# fig, ax = plt.subplots()
# ax.set_title("BCC", size=15)
# ax.scatter(k["k"], k["delta_k"], facecolors="none", edgecolors="k", marker="*", s=95)
# ax.set_xlabel(r"K (nm$^{-1}$)", size=15)
# ax.set_ylabel(r"$\Delta$K (nm$^{-1}$)", size=15)
# ax.tick_params(axis="both", labelsize=12)
# ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)

# for i in range(len(planes)):
#     ax.annotate("{" + planes[i] + "}", (1.01*k["k"][i], 1.01*k["delta_k"][i]), fontsize=12)
# ax.set_xlim(0, .12)
# ax.set_ylim(0., .12)

# fig.tight_layout()







# fig1, ax1 = plt.subplots()
# ax1.set_title("BCC", size=15)
# ax1.scatter(h_squared, (delta_k.T[100] - 0.00000)**2 / k.T[100]**2, facecolors="none", edgecolors="k", marker="*", s=95)
# ax1.set_xlabel(r"H$^2$", size=15)
# ax1.set_ylabel(r"$\frac{\left(\Delta K - \alpha \right)^2}{k^2}$", size=15)
# ax1.tick_params(axis="both", labelsize=12)
# ax1.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)

# print(delta_k.T)

# print()
# print(delta_k)
# print()
# print()
# for i in range(len(planes)):
#     ax1.annotate("{" + planes[i] + "}", (1.01*h_squared[i], 1.01*y_values["y"][i]), fontsize=12)
    
# ax1.annotate(f"R²: {y_values['r_sq']:.4f}", (.25, y_values["y"].max()), fontsize=12)
# ax1.annotate(rf"$\alpha$: {y_values['alpha']}", (.25, y_values["y"].max()*.975), fontsize=12)

# # ax.set_xlim(0, .12)
# # ax.set_ylim(0., .2)

# x1 = np.arange(-.02, .35, .01)
# ax1.plot(x1, y_values["coef"] * x1 + y_values["intercept"])
# ax1.annotate(rf"q: {- float(y_values['coef']) / float(y_values['intercept']):.3f}", (.25, y_values["y"].max()*.95), fontsize=12)


# fig1.tight_layout()

# plt.show()