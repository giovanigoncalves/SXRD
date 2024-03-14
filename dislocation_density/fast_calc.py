import numpy as np 
import pandas as pd


 
c11 = 276.6
c12 = 145.8
c44 = 97.58

a_q_screw = 6.0725
b_q_screw = 0.4338
c_q_screw = 0.0415
d_q_screw = -3.5021

a_q_edge = 8.8331
b_q_edge = 0.8241
c_q_edge = 0.1078
d_q_edge = -7.0570


A = (2 * c44) / (c11 - c12)

q_screw = a_q_screw * (1 - np.exp(-A / b_q_screw)) + c_q_screw * A + d_q_screw

q_edge = a_q_edge * (1 - np.exp(-A / b_q_edge)) + c_q_edge * A + d_q_edge

# print("\n\n")
# print(f"A: {A:.2f}")
# print(f"Screw: {q_screw:.2f}")
# print(f"Edge: {q_edge:.2f}")
# print("\n\n")

if 1 < 2 & 3 < 4:
    print("ok")