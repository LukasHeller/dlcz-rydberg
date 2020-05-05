import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

p = {"bsplit": 0.47,        # Beam splitter ratio
     "eta_exc": 0.45,       # Excelitas quantum efficiency
     "eta_count": 0.65,     # Laseroptics Count quantum efficiency
     "eta_D_int": 0.35,     # DLCZ intrinsic retrieval efficiency
     "eta_D_cou_r": 0.85,   # DLCZ read coupling efficiency
     "eta_D_tr_r": 0.9,     # DLCZ read transmission efficiency
     "eta_D_cou_w": 0.75,   # DLCZ write coupling efficiency
     "eta_D_tr_w": 0.9,     # DLCZ transmission efficiency
     "eta_D_det_w": 0.45,   # DLCZ write detection efficiency
     "D_R": 1/10,           # Rydberg duty cicle
     "D_D": 1.5/(9+1.5),    # DLCZ duty cicle
     "R_D": 1/(2*10**-6),   # Write trial repetition rate
     "eta": 0.9,
     "p_R": 0.04}           # Photon overlap factor

with open('p.json', 'w') as fp:
    json.dump(p, fp)

W = np.linspace(1,100,100)      # Detected DLCZ write photon rate per second

#%% Simulation function definition

def sim(p,W):
    d = pd.DataFrame(index = W)
    
    # Multiplicative detection efficiency
    eta_det = np.sqrt(p["bsplit"]*(1-p["bsplit"])*p["eta_exc"]*p["eta_count"])
    
    # Probability to find a heralded DLCZ read photon in front of the beam splitter 
    p_D = p["eta_D_int"]*p["eta_D_cou_r"]*p["eta_D_tr_r"]
    
    # Probability to find a "heralded" Rydberg photon in front of the bea, splitter
    #p_R = p_R_gen*eta_R_cou_r*eta_R_tr_r*eta_R_AOM
    p_R = p["p_R"]
    
    # DLCZ Write photon detection probability per trial
    p_w_det = W/(p["D_R"]*p["D_D"]*p["R_D"])
    
    # DLCZ write photon emission probability per trial
    p_D_ext = p_w_det/(p["eta_D_cou_w"]*p["eta_D_tr_w"]*p["eta_D_det_w"])
    
    # DLCZ cross correlation function g_wr
    g_D_wr = 1 +(p["eta_D_int"]*(1-p_D_ext))/(p_D_ext*p["eta_D_int"] + p_D_ext*(1-p["eta_D_int"])*0.5)
    
    # DLCZ autocorrelation function g_rr,w, based on the crosscorrelation
    #g_D = (4*W)/(p["eta_D_cou_w"]*p["eta_D_tr_w"]*p["eta_D_det_w"]*p["D_R"]*p["D_D"]*p["R_D"])
    g_D = 4/g_D_wr

    # Rydberg autocorrelation function g_rr
    g_R = 0.25
    
    # Probability of coincidence click per trial, indistinguishable case
    p_12 = p_R**2*eta_det**2*g_R + p_D**2*eta_det**2*g_D + 2*eta_det**2*p_R*p_D*(1-p["eta"])
    
    # Probability of coincidence click per trial, distinguishable case
    p_12_0 = p_R**2*eta_det**2*g_R + p_D**2*eta_det**2*g_D + 2*eta_det**2*p_R*p_D
    
    # Visibility
    V = 1-p_12/p_12_0
    
    # Mysterious parameter of coincidence count rate in distinghuishable case.
    # --- I GUESS THIS IS VERY WRONG ??? ---
    X = d.index*p_12_0
    
    # Append to dataframe
    d["g_D"] = g_D.copy()
    d["V"] = V.copy()
    d["p_w_det"] = p_w_det.copy()
    d["X"] = X.copy()
    d["p_12_0"] = p_12_0.copy()
    d["p_12"] = p_12.copy()
    
    return d

d = sim(p,W)
#%% Some simulations

fig, ax = plt.subplots(1,4, num = 5, figsize = (12,4))
p_c = p.copy()

fig1, ax1 = plt.subplots()
for p_R in np.linspace(0.01,1,100):
     p_c["p_R"] = p_R
     
     d = sim(p_c, np.array([10]))
     ax1.scatter(p_R,d.V)
#%%
# for bsplit in [0.3,0.40, 0.5]:
#     p_c["bsplit"] = bsplit
#     label = p_c["bsplit"]

# for eta_D_int in [0.25,0.35, 0.45]:
#     p_c["eta_D_int"] = eta_D_int
#     label = p_c["eta_D_int"]

#for D_R in [0.25,0.35, 0.45]:
#    p_c["D_R"] = D_R
#    label = p_c["D_R"]

#for eta_D_cou_r in [0.4, 0.6, 0.8]:
#    p_c["eta_D_cou_r"] = eta_D_cou_r
#    label = p_c["eta_D_cou_r"]
#
for p_R in [0.04, 0.06, 0.08]:
    p_c["p_R"] = p_R
    label = p_c["p_R"]
    
    d = sim(p_c,W)
#    print(d["X"])
#    print(d["V"])
    
#    ax[0].plot(X*3600,d["V"], label = label)
#    ax[0].set_xlabel("Coincidences distinghuishable case (1/h)")
#    ax[0].set_ylabel("Visibility")
#    
#    ax[3].plot(X*3600,d["p_12_0"], label = label)
#    ax[3].set_xlabel("Coincidences distinghuishable case (1/h)")
#    ax[3].set_ylabel("p_12_0")
#    
#    ax[1].plot(X*3600,d["g_D"], label = label)
#    ax[1].set_xlabel("Coincidences distinghuishable case (1/h)")
#    ax[1].set_ylabel(r"$g_{r,r}^{(2)}$ DLCZ")
#    
#    ax[2].plot(X*3600,d["p_w_det"]*100, label = label)
#    ax[2].set_xlabel("Coincidences distinghuishable case (1/h)")
#    ax[2].set_ylabel("p_w_det (%)")
    
    ax[0].plot(d["X"]*3600,d["V"], label = label)
    ax[0].set_xlabel("Coincidences distinghuishable case (1/h)")
    ax[0].set_ylabel("Visibility")
    
    ax[3].plot(d["X"]*3600,d["p_12_0"], label = label)
    ax[3].set_xlabel("Coincidences distinghuishable case (1/h)")
    ax[3].set_ylabel("p_12_0")
    
    ax[1].plot(d["X"]*3600,d["g_D"], label = label)
    ax[1].set_xlabel("Coincidences distinghuishable case (1/h)")
    ax[1].set_ylabel(r"$g_{r,r}^{(2)}$ DLCZ")
    
    ax[2].plot(d["X"]*3600,d["p_w_det"]*100, label = label)
    ax[2].set_xlabel("Coincidences distinghuishable case (1/h)")
    ax[2].set_ylabel("p_w_det (%)")
    
ax[0].legend(title = "p_R")
plt.tight_layout()
plt.show()

plt.savefig("sim.png", dpi = 300, bbox_inches = "tight")

 


