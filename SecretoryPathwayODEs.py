import numpy as np
from scipy.integrate import odeint
from math import cos,pi,exp
import matplotlib.pyplot as plt

# Define species concentrations
#C_ER                 # μM        ER-located collagen
#C_H                  # μM        Collagen-HSP47 complex in ERGIC
#C_G                  # μM        Golgi-located collagen
#C_PG                 # μM        Post-Golgi compartments-located collagen
#C_PM                 # μM        Plasma membrane-located collagen
#C_E                  # μM        Extracellular collagen
#H_ER                 # μM        ER-located HSP47
#H_G                  # μM        Golgi-located HSP47
#K                    # μM        Pka
#S                    # μM        SEC61A2
#T                    # μM        TANGO1
#P                    # μM        PDE4D
#V                    # μM        VPS33B
#M                    # μM        MMP14
#D                    # μM        CTSK

# Define mass action kinetics parameters
μ_typ    = 0.0342    # μM        Median typical concentration
k_1      = 0.00822   # μM        SEC61A2 mean concentration
k_2      = 0.0105    # μM        TANGO1 mean concentration
k_3      = 0.0342    # μM        PDE4D mean concentration
k_6      = 0.00335   # μM        VPS33B mean concentration
k_7      = 0.00457   # μM        MMP14 mean concentration
k_15     = 0.0342    # μM        CTSK mean concentration
t_bar    = 24.0      # h         Circadian period
k_4a     = 4.84      # h⁻¹       HSP47 protein synthesis rate
k_4b     = 0.00413   # h⁻¹       ER HSP47 degradation rate
k_4c     = 0.00413   # h⁻¹       Golgi HSP47 degradation rate
k_5a     = 1.89      # h⁻¹       PKA protein synthesis rate
k_5b     = 0.0069    # h⁻¹       PKA protein degradation rate
k_8      = 0.0693    # h⁻¹       Col-I synthesis rate
k_10     = 3600.0    # h⁻¹       ERGIC to Golgi Col-I transition rate
k_12     = 0.0417    # h⁻¹       Post-golgi compartments to plasma membrane Col-I transition rate
k_5c     = 29.2      # μM⁻¹      Rate constant for repression of PKA by PDE4D
k_minus9 = 1.22      # μM⁻¹h⁻¹   Col-I-HSP47 dissociation rate
k_11     = 1.22      # μM⁻¹h⁻¹   Golgi to post-golgi compartments Col-I transition rate
k_13     = 12.2      # μM⁻¹h⁻¹   Plasma membrane to extracellular Col-I transition rate
k_14     = 1.22      # μM⁻¹h⁻¹   CTSK-dependent extracellular Col-I degradation rate
k_16     = 1.22      # μM⁻¹h⁻¹   Golgi to ER HSP47 transition rate
k_plus9  = 0.356     # μM⁻²h⁻¹   Col-I-HSP47 complex association rate
G_bar    = 2.0       # Unitless  Mean pulse function value determined by mean value theorem
#t_0                  # h         Phase


def F(t,t_0):
    output = 0.5*cos(2.0*pi*(t-t_0)/t_bar)+1.0
    return output


#cos_e = cos(exp(1))
#cos_minus_e = cos(exp(-1))
def G(t,t_0):
    timedependentcomponent = cos(exp(cos(2.0*pi*(t-t_0)/t_bar)))
    numerator = cos(exp(-1.0)) - timedependentcomponent
    denominator = cos(exp(-1.0)) - cos(exp(1))
    output = ((numerator/denominator)**4)/G_bar
    return output

def SecretoryPathwayODEs(y,t):
    C_ER = y[0]
    C_H  = y[1]
    C_G  = y[2]
    C_PG = y[3]
    C_PM = y[4]
    C_E  = y[5]
    H_ER = y[6]
    H_G  = y[7]
    K    = y[8]

    S = k_1*F(t,3.0)
    T = k_2*G(t,9.0)
    P = k_3*F(t,19.0)
    V = k_6*F(t,19.0)
    M = k_7*G(t,3.0)
    D = k_15*F(t,9.0)


    dC_ERdt     = k_8*S - k_plus9*C_ER*T*H_ER + k_minus9*C_H*T
    dC_Hdt      = k_plus9*C_ER*T*H_ER - k_minus9*C_H*T - k_10*C_H
    dC_Gdt      = k_10*C_H - k_11*C_G*V
    dC_PGdt     = k_11*C_G*V - k_12*C_PG
    dC_PMdt     = k_12*C_PG - k_13*C_PM*M
    dC_Edt      = k_13*C_PM*M - k_14*D*C_E
    dH_ERdt     = k_4a*S - k_4b*H_ER - k_plus9*C_ER*T*H_ER + k_minus9*C_H*T + k_16*H_G*K
    dH_Gdt      = k_10*C_H - k_16*H_G*K - k_4c*H_G
    dKdt        = k_5a*S/(1+k_5c*P) - k_5b*K

    x = np.array([dC_ERdt,dC_Hdt,dC_Gdt,dC_PGdt,dC_PMdt,dC_Edt,dH_ERdt,dH_Gdt,dKdt])

    return x


# Initial conditions
ndays = 160
y0 = np.zeros(9)
t = np.linspace(0,24*ndays,1000*ndays)

output = odeint(SecretoryPathwayODEs,y0,t[1:])

names = ["C_ER","C_H ","C_G ","C_PG","C_PM","C_E ","H_ER","H_G ","K   "]

fig,ax = plt.subplots()

for i in range(np.shape(y0)[0]):
    ax.plot(t[1:],output[:,i],label=names[i])

ax.set_ylim([0,10])
ax.set_xlim([0,3600])
ax.legend()
plt.savefig("/Users/christopher/Desktop/test.png",dpi=500,bbox_inches='tight',padding_inches=0)
plt.show()
