import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from sklearn.linear_model import LinearRegression
import os
from pprint import pprint
import json
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

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
df = pd.read_excel("Juan_diffraction.xlsx", skiprows=1, names=names)
df.replace("NP", np.nan, inplace=True)


# Chosing one unique diffraction spectra to calc dislocation density (FWHM and angle info)
SPECTRAS = np.arange(1, 15, 1) # Number of the diffraction data
INSTRUMENTAL_BROADENING = 0.0160558308639235

fwhm = {}
angle = {}

struc = "FCC"
if struc == "BCC":
    for i in SPECTRAS:
        fwhm[i] = [df.M110_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.M200_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.M211_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.M220_fwhm[i] - INSTRUMENTAL_BROADENING]
        fwhm[i] = np.array([float(n) for n in fwhm[i]])
        fwhm[i] = np.array([np.deg2rad(n) for n in fwhm[i]]) # Transforming the fwhm from degree to radians
        
        angle[i] = [df.M110_angle[i], df.M200_angle[i], df.M211_angle[i], df.M220_angle[i]]
        angle[i] = np.array([np.deg2rad(n/2) for n in angle[i]]) # Transforming the angles from degree to radians

elif struc == "FCC":
    for i in SPECTRAS:
        fwhm[i] = [df.A111_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.A200_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.A220_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.A311_fwhm[i] - INSTRUMENTAL_BROADENING, 
                   df.A222_fwhm[i] - INSTRUMENTAL_BROADENING]
        fwhm[i] = np.array([float(n) for n in fwhm[i]])
        fwhm[i] = np.array([np.deg2rad(n) for n in fwhm[i]]) # Transforming the fwhm from degree to radians
        
        angle[i] = [df.A111_angle[i], df.A200_angle[i], df.A220_angle[i], df.A311_angle[i], df.A222_angle[i]]
        angle[i] = np.array([np.deg2rad(n/2) for n in angle[i]]) # Transforming the angles from degree to radians

        
# Here we start the play

class DislocationDensityMWH:
    
    # Initial data (wavelength and phase planes)
    def __init__(self, structure):
        self.wavelength = 1.4235 # in nanometers
        
        self.structure = structure
        if self.structure == "BCC":
            self.planes = ["110", "200", "211", "220"] # for BCC structure
        elif self.structure == "FCC":
            self.planes = ["111", "200", "220", "311", "222"] # for FCC structure
            self.Wg = [0.43, 1.0, 0.71, 0.45, 0.43]
        self.a_bcc = 2.87
        self.a_fcc = 3.59
        self.alpha_range = np.arange(0, 1e-4, 1e-9)
        self.beta_prime_range = np.arange(0, 1e-3, 1e-7)
        
    # Calculating the K values for the different planes
    def k_values(self, angle, fwhm):
        """
            Calculates the K and Delta_K values. 
        """
        k = {}
        delta_k = {}
        for key, value in angle.items():     
            k[key] = 2 * np.sin(value) / self.wavelength # the angles need to be in radians
            
        for i in SPECTRAS:
            delta_k[i] = (2 * np.cos(angle[i]) * fwhm[i]) / self.wavelength
            
        # print(k)
        # print(delta_k)
        return dict(k=k, delta_k=delta_k)
    

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
        print(h_values)
        return h_values


    def alpha_value(self, angle, fwhm):
        k = self.k_values(angle, fwhm)
        h_squared = self.H_squared().values() # x values for the plot
        h_squared = [i for i in h_squared]
        x = [i for i in h_squared]
        x = np.array(x).reshape((-1, 1))
        new_x = np.arange(0, 1, .01)

        best_alpha = {}
        best_r_sq = {}
        coef = {}
        intercept = {}
        q = {}
        

        p = 0
        for key, value in k["delta_k"].items():
            best_alpha[key] = 0
            best_r_sq[key] = 0
            coef[key] = 0
            intercept[key] = 0
            for num in self.alpha_range:
                y =  (k["delta_k"][key]**2 - num) / k["k"][key]**2
                
                model = LinearRegression().fit(x, y)
                r_sq = model.score(x, y)
                if (r_sq > best_r_sq[key]) & (r_sq < 1):
                    best_r_sq[key] = r_sq
                    best_alpha[key] = num
                    coef[key] = model.coef_
                    intercept[key] = model.intercept_
                p += 1
                print(f"{p*100/(len(self.alpha_range)*len(k['delta_k'].keys())):.2f} %")
            q[key] = float(( -1 * coef[key] ) / intercept[key])            
            
            fig, ax = plt.subplots()
            y_value = np.array((k["delta_k"][key]**2 - best_alpha[key]) / k["k"][key]**2 )
            y_err = 2 * (1 - np.array(best_r_sq[key])) * y_value
            

            ax.errorbar(x=h_squared, 
                        y=y_value, 
                        yerr=y_err,
                        c="k",
                        capsize=7,
                        linestyle="",
                        marker="o",
                        mfc="white")
    

        return dict(alpha=best_alpha, r_sq=best_r_sq, coef=coef, intercept=intercept, q=q)    
    
    
    def q_values(self, angle, fwhm):
        
        alpha_values = self.alpha_value(angle, fwhm)
        alpha = alpha_values["alpha"]
        q = alpha_values["q"]
        r_sq = alpha_values["r_sq"]
        intercept = alpha_values["intercept"]
        
        if self.structure.upper() == "BCC":
            df_screw = pd.read_csv(f"abcd_parameters_for_q_screw_dislocation_bcc.txt", delimiter=";")
            df_edge = pd.read_csv(f"abcd_parameters_for_q_edge_dislocation_bcc.txt", delimiter=";")
            # c11, c12 and c44 were chosen from "Elastic Properties of Fe–C and Fe–N Martensites"
            # Maaouia SOUISSI and Hiroshi NUMAKURA
            # ISIJ International, Vol. 55 (2015), No. 7, pp. 1512–1521
            
            c11 = 276.6
            c12 = 145.8
            c44 = 97.58
            
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

   
        elif self.structure.upper() == "FCC":
            df_edge = pd.read_csv(f"abcd_parameters_for_q_edge_dislocation_fcc.txt", delimiter=";")
            # c11, c12 and c44 were chosen from F. HajyAkbary et al. (2015) Materials Science & Engineering A 639
            c11 = 163
            c12 = 141
            c44 = 79.6
            
            a_q_screw = 5.4252
            b_q_screw = 0.7196
            c_q_screw = 0.0690
            d_q_screw = -3.1970
            
            if abs(c12/c44 - 0.5) < abs(c12/c44 - 1) < abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_q_edge = df_edge["c12/c44=0.5"][0]
                b_q_edge = df_edge["c12/c44=0.5"][1]
                c_q_edge = df_edge["c12/c44=0.5"][2]
                d_q_edge = df_edge["c12/c44=0.5"][3]

            elif abs(c12/c44 - 1) < abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_q_edge = df_edge["c12/c44=1"][0]
                b_q_edge = df_edge["c12/c44=1"][1]
                c_q_edge = df_edge["c12/c44=1"][2]
                d_q_edge = df_edge["c12/c44=1"][3]

            elif abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_q_edge = df_edge["c12/c44=2"][0]
                b_q_edge = df_edge["c12/c44=2"][1]
                c_q_edge = df_edge["c12/c44=2"][2]
                d_q_edge = df_edge["c12/c44=2"][3]
                
            else:
                a_q_edge = df_edge["c12/c44=3"][0]
                b_q_edge = df_edge["c12/c44=3"][1]
                c_q_edge = df_edge["c12/c44=3"][2]
                d_q_edge = df_edge["c12/c44=3"][3]
            
        else:
            print("Invalid structure! Chose BCC or FCC!")
                   

        A = (2 * c44) / (c11 - c12) # Elastic anisotropy
        q_screw = a_q_screw * (1 - np.exp(-A / b_q_screw)) + c_q_screw * A + d_q_screw # constant value for 100% screw dislocations
        q_edge = a_q_edge * (1 - np.exp(-A / b_q_edge)) + c_q_edge * A + d_q_edge # constant value for 100% edge dislocations

        
        f_edge = {}
        f_screw = {}
        for key in q:
            if q[key] < q_edge:
                f_edge[key] = 1
                f_screw[key] = 0
            elif q[key] > q_screw:
                f_edge[key] = 0
                f_screw[key] = 1
            else:
                f_edge[key] = ( q_screw - q[key] ) / ( q_screw - q_edge )  
                f_screw[key] = 1 - f_edge[key]

            print("\n\n")
            
            
            print(f"q_screw: {q_screw}")
            print(f"q_edge: {q_edge}")
            print(f"q_exp: {q[key]}")
            
            
            print("\n\n")
            print(f"f_q_screw for {key}: {100 * float(f_screw[key]):.1f} %")
            print(f"f_q_edge for {key}: {100 * float(f_edge[key]):.1f} %")
            print("\n\n")
            
        return dict(f_edge=f_edge, f_screw=f_screw, 
                    q_screw=q_screw, q_edge=q_edge, q=q, 
                    A=A, r_sq=r_sq, intercept=intercept,
                    alpha=alpha)

    
    def ch00(self, angle, fwhm):
        
        q = self.q_values(angle, fwhm)
        alpha = q["alpha"]
    
        if self.structure.upper() == "BCC":
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
                a_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][3])

            elif abs(c12/c44 - 1) < abs(c12/c44 - 2):
                a_ch00_edge = float(df_ch00_edge["c12/c44=1"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=1"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=1"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=1"][3])

            else:
                a_ch00_edge = float(df_ch00_edge["c12/c44=2"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=2"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=2"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=2"][3])

        
        elif self.structure.upper() == "FCC":
            df_ch00_edge = pd.read_csv("abcd_parameters_for_ch00_edge_dislocation_fcc.txt", delimiter=";")
            # Ref. G. Ghosh and G. B. Olson. Acta Materialia 50 (2002) 2655-2675
            c11 = 540
            c12 = 120
            c44 = 260
            
            # Os valores praticamente independem do c12/c44 para screw dislocations in FCC structure
            a_ch00_screw = 0.174
            b_ch00_screw = 1.9522
            c_ch00_screw = 0.0293
            d_ch00_screw = 0.0662

            if abs(c12/c44 - 0.5) < abs(c12/c44 - 1) < abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=0.5"][3])

            elif abs(c12/c44 - 1) < abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_ch00_edge = float(df_ch00_edge["c12/c44=1"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=1"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=1"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=1"][3])
                
            elif abs(c12/c44 - 2) < abs(c12/c44 - 3):
                a_ch00_edge = float(df_ch00_edge["c12/c44=2"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=2"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=2"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=2"][3])

            else:
                a_ch00_edge = float(df_ch00_edge["c12/c44=3"][0])
                b_ch00_edge = float(df_ch00_edge["c12/c44=3"][1])
                c_ch00_edge = float(df_ch00_edge["c12/c44=3"][2])
                d_ch00_edge = float(df_ch00_edge["c12/c44=3"][3])
        
        ch00_edge = a_ch00_edge * ( 1 - np.exp( - float(q["A"]) / b_ch00_edge) ) + c_ch00_edge * float(q["A"]) + d_ch00_edge
        ch00_screw = a_ch00_screw * ( 1 - np.exp( - float(q["A"]) / b_ch00_screw) ) + c_ch00_screw * float(q["A"]) + d_ch00_screw
        
        ch00 = {key: float(q["f_edge"][key]) * ch00_edge + float(q["f_screw"][key]) * ch00_screw for key, value in q["f_edge"].items()}
 
        return dict(ch00_edge=ch00_edge, ch00_screw=ch00_screw, ch00=ch00, q=q, alpha=alpha)
        
        
    def c(self, angle, fwhm):
        ch00 = self.ch00(angle, fwhm)
        alpha = ch00["alpha"]
        h_values = self.H_squared()
        c = {}
        for key, value in ch00["ch00"].items():
            l = []
            for j in self.planes:
                l.append(ch00["ch00"][key] * ( 1 - ch00["q"]["q"][key] * h_values[j]))
            c[key] = l
        pprint(c)
        return c, ch00["q"]["r_sq"], ch00["q"]["intercept"], ch00["q"]["f_edge"], ch00["q"]["f_screw"], alpha
     
    
    def b(self):

        b = {}
        if self.structure.upper() == "BCC":
            a = self.a_bcc
            b["111"] = (a / 2) * np.sqrt(1**2 + 1**2 + 1**2)
            print(self.structure)
            return b["111"] # According to ref. T. Ungár et al. (1999) pag. 993 {110}<111>
        elif self.structure.upper() == "FCC":
            a = self.a_fcc
            b["110"] = (a / 2) * np.sqrt(1**2 + 1**2 + 0**2)
            print(self.structure)
            return b["110"] # According to ref. T. Ungár et al. (1999) pag. 993 {111}<110>

        
    def dislocation_density(self, angle, fwhm, M=1.4):
        '''
            M value according to ref: F. HajyAkbary et al. (2015)
        '''
        
        b = self.b()
        c, r_sq_1, intercept_grain, f_edge, f_screw, alpha = self.c(angle, fwhm)
        k = self.k_values(angle, fwhm)

        if self.structure.upper() == "BCC":
            coef = {}
            intpt = {}
            r_sq_2 = {}
            fig, ax = plt.subplots()
            for key, value in k["k"].items():
            
                x = np.array(k["k"][key]) * np.sqrt(np.array(c[key]))
                # print(x)    
                x_ = np.array(x).reshape((-1, 1))
                # print(x_)
                new_x = np.linspace(x.min(), x.max(), 1000)
                
                
                y =  k["delta_k"][key]
                # print(y)
                try:
                    model = LinearRegression().fit(x_, y) 
                    r_sq_2[key] = model.score(x_, y)
                    intpt[key] = model.intercept_
                    coef[key] = model.coef_
                    
                    fig, ax = plt.subplots()
                    y_err = 2 * (1 - np.array(r_sq_2[key])) * y
                    ax.errorbar(x=x, 
                            y=y, 
                            yerr=y_err,
                            c="k",
                            capsize=7,
                            linestyle="",
                            marker="o",
                            mfc="white")
                    # ax.plot(x, y, c="k", marker="o", mfc="white", ls="")
                    ax.plot(new_x, intpt[key] + coef[key]*new_x, c="k")
                    ax.set_xlabel(r"K$\overline{C}$$^{1/2}$ (nm$^{-1}$)", size=15)
                    ax.set_ylabel(r"$\Delta K$ (nm$^{-1}$)", size=15)
                    ax.tick_params(axis="both", labelsize=12)
                    # ax.set_ylim(0, )
                    # ax.annotate(f"Measure: {key}", (0, .01))
                    ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)
                    ax.minorticks_on()
                    fig.tight_layout()
                
                    try:
                        os.mkdir(f"dislocation_density_MWH_method_plot_{self.structure}")
                    except FileExistsError:
                        pass
                    
                    plt.savefig(f"dislocation_density_MWH_method_plot_{self.structure}/Measure_{50 + 100 * (key)}.png")
                except ValueError:
                    pass
                       
                
            def rho_f(coef, b, M):
                return 1e20 * 2 * coef**2 / (b**2 * M**2 * np.pi) #in 1/m²
        
            def std_value(value, r_sq_1, r_sq_2):
                return 2 * value * (1 - (r_sq_1 * r_sq_2))
            
            rho = {key: rho_f(coef[key], b, M) for key, value in coef.items()}
            
            rho_value = {key: {"rho": rho[key], 
                            "std": std_value(rho[key], r_sq_1[key], r_sq_2[key]),
                            "D": 0.9 / intercept_grain[key],
                            "std_grain": std_value(intercept_grain[key], r_sq_1[key], r_sq_2[key]),
                            "f_edge": f_edge[key],
                            "f_screw": f_screw[key],
                            "alpha": alpha[key]
                            } for key, value in rho.items()}
        
            names = [i for i in rho_value.keys()]
            values = [[rho_value[col]["D"], 
                    rho_value[col]["std_grain"], 
                    float(rho_value[col]["rho"]), 
                    float(rho_value[col]["std"]),
                    float(rho_value[col]["f_edge"]),
                    float(rho_value[col]["f_screw"]),
                    float(rho_value[col]["alpha"])
                    ] for col in names]
            files = [i+1 for i in names]
            # print(values)
            cols = ["D", "std_grain", "rho", "std", "f_edge", "f_screw", "alpha"]
        
        
        elif self.structure.upper() == "FCC":
            Wg = np.array([.43, 1.0, .71, .45, .43]) # for 111, 200, 220, 311, and 222 austenite planes, respectively

            best_beta = {}
            r_sq_2 = {}
            a_coef = {}
            b_coef = {}
            c_coef = {}
            
            def func(x, a, b, c):
                return a*x**2 + b*x + c
            
            p = 0
            for key, value in k["delta_k"].items():
                x = np.array(k["k"][key]) * np.sqrt(np.array(c[key]))
                x_ = np.array(x).reshape((-1, 1))
                new_x = np.linspace(x.min(), x.max(), 1000)
                
                best_beta[key] = 0
                r_sq_2[key] = 0
                a_coef[key] = 0
                b_coef[key] = 0
                c_coef[key] = 0
                for num in self.beta_prime_range:
                    y =  np.array(k["delta_k"][key] - num * Wg)
                    # print(y)
                    popt, pcov = curve_fit(func, x, y)
                    
                    y_pred = func(x, *popt)
                    r_squared = r2_score(y, y_pred)
                    print(r_squared)
                    if r_squared > r_sq_2[key]:
                        r_sq_2[key] = r_squared
                        best_beta[key] = num
                        a_coef[key] = popt[0]
                        b_coef[key] = popt[1]
                        c_coef[key] = popt[2]
                    p += 1
                    print(f"{p*100/(len(self.beta_prime_range)*len(k['delta_k'].keys())):.2f} %")
                print(best_beta)
                print(r_sq_2)
                
                fig, ax = plt.subplots()
                y_value = np.array(k["delta_k"][key] - best_beta[key] * Wg)
                y_err = 2 * (1 - np.array(r_sq_2[key])) * y_value

                ax.errorbar(x, y_value, yerr=y_err, c="k", marker="o", mfc="white", ls="", capsize=5)
                ax.plot(new_x, func(new_x, a_coef[key], b_coef[key], c_coef[key]), c="k")
                ax.set_xlabel(r"K$\overline{C}$$^{1/2}$ (nm$^{-1}$)", size=15)
                ax.set_ylabel(r"$\Delta K - \beta 'W_{hkl}$ (nm$^{-1}$)", size=15)
                ax.tick_params(axis="both", labelsize=12)
                ax.set_ylim(0, )
                # # ax.annotate(f"Measure: {key}", (0, .01))
                ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0), useMathText=True)
                ax.minorticks_on()
                fig.tight_layout()
                # plt.show()
                
                try:
                    os.mkdir(f"dislocation_density_MWH_method_plot_{self.structure}")
                except FileExistsError:
                    pass
                
                plt.savefig(f"dislocation_density_MWH_method_plot_{self.structure}/Measure_{50 + 100 * (key)}_{self.structure}.png")
        
            def rho_f(b_coef, b, M):
                return 1e20 * 2 * b_coef**2 / (b**2 * M**2 * np.pi) #in 1/m²
            
            def std_value(value, r_sq_1, r_sq_2):
                return 2 * value * (1 - (r_sq_1 * r_sq_2))
            
            rho = {key: rho_f(b_coef[key], b, M) for key, value in b_coef.items()}
            
            rho_value = {key: {"rho": rho[key], 
                            "std": std_value(rho[key], r_sq_1[key], r_sq_2[key]),
                            "D": 0.9 / c_coef[key],
                            "std_grain": std_value(c_coef[key], r_sq_1[key], r_sq_2[key]),
                            "f_edge": f_edge[key],
                            "f_screw": f_screw[key],
                            "beta_prime": best_beta[key],
                            "alpha": alpha[key]
                            } for key, value in rho.items()}
        
            names = [i for i in rho_value.keys()]
            values = [[rho_value[col]["D"], 
                    rho_value[col]["std_grain"], 
                    float(rho_value[col]["rho"]), 
                    float(rho_value[col]["std"]),
                    float(rho_value[col]["f_edge"]),
                    float(rho_value[col]["f_screw"]),
                    float(rho_value[col]["beta_prime"]),
                    float(rho_value[col]["alpha"])
                    ] for col in names]
            files = [i+1 for i in names]
            # print(values)
            cols = ["D", "std_grain", "rho", "std", "f_edge", "f_screw", "beta_prime", "alpha"]
      
         
            
        # Thus...
        df = pd.DataFrame(values, columns=cols)
        df["names"] = files
        df.set_index("names", inplace=True)
        print(df)

        df.to_csv(f"dislocation_density_data_{self.structure}.txt", sep=",")

        pprint(rho_value)
        return rho_value 


a = DislocationDensityMWH("FCC")
# a.structure = "FCC"
# k = a.k_values(angle, fwhm)
# k = k.to_dict()
# h_squared = a.H_squared()
# h_squared = [value for value in h_squared.values()]
# alpha = a.alpha_value(angle, fwhm)
# pprint(alpha)
# q = alpha["q"]
# q_values = a.q_values(angle, fwhm)
# ch00 = a.ch00(angle, fwhm)
# c = a.c(angle, fwhm)

dis = a.dislocation_density(angle, fwhm)



# pprint(dis)
