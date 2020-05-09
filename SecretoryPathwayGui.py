import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button
from scipy.integrate import odeint
from math import cos,pi,exp






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

cos_e = cos(exp(1.0))
cos_minus_e = cos(exp(-1.0))
def G(t,t_0):
    timedependentcomponent = cos(exp(cos(2.0*pi*(t-t_0)/t_bar)))
    numerator = cos_minus_e - timedependentcomponent
    denominator = cos_minus_e - cos_e
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






initA = 1.0
initB = 5.0

# Set up plot canvas
fig,ax=plt.subplots()
plt.subplots_adjust(left=0.1,bottom=0.35)
#ax.set_xlim([0,96])
ax.set_ylim([0,0.2])

#Create slider to vary A value
axSlider1 = plt.axes([0.1,0.2,0.8,0.05])
axSlider1.set_xticks([0.0,100.0])
slider1 = Slider(axSlider1,'a',valinit=initA,valmin=0,valmax=5)
axSlider1.xaxis.set_visible(True)
axSlider1.set_xticks([0,5])

#Create slider to vary B value
axSlider2 = plt.axes([0.1,0.1,0.8,0.05])
axSlider2.set_xticks([0.0,100.0])
slider2 = Slider(axSlider2,'b',valinit=initB,valmin=0,valmax=10,color='orange')
axSlider2.xaxis.set_visible(True)
axSlider2.set_xticks([0,10])

# Global variable to switch dynamic plotting on and off
onOff = False
time = 0
dt = 0.1
# Initialise a flat line
x = np.arange(0,96,dt)
y = np.zeros([np.shape(x)[0],9])
#Init_vals = np.zeros(9)
output = np.zeros([20,9])
line0, = ax.plot(x,y[:,0])
line1, = ax.plot(x,y[:,1])
line2, = ax.plot(x,y[:,2])
line3, = ax.plot(x,y[:,3])
line4, = ax.plot(x,y[:,4])
line5, = ax.plot(x,y[:,5])
line6, = ax.plot(x,y[:,6])
line7, = ax.plot(x,y[:,7])
line8, = ax.plot(x,y[:,8])

# Function to continually update the plot so long as onOff=True
def updateLine(val):
    global onOff
    global time
    global y
    global output
    # Set onOff = True if currently False
    if onOff:
        pass
    else:
        onOff = True
    # Extract A and B values from sliders
    a_val = slider1.val
    b_val = slider2.val
    while onOff:
        trange = np.linspace(time,time+20*dt,num=20,endpoint=False)
        output = odeint(SecretoryPathwayODEs,y[-1,:],trange)
        time = time+20*dt
        # Shift most of the values left by 20 points
        y[:-20]=y[20:]
        # Fill the last 20 points at the right of the plot with Gaussian noise adjusted with A and B
        y[-20:,:]=output[:,:]
        # Redraw line
        line0.set_ydata(y[:,0])
        line1.set_ydata(y[:,1])
        line2.set_ydata(y[:,2])
        line3.set_ydata(y[:,3])
        line4.set_ydata(y[:,4])
        line5.set_ydata(y[:,5])
        line6.set_ydata(y[:,6])
        line7.set_ydata(y[:,7])
        line8.set_ydata(y[:,8])
        plt.draw()
        # Stop briefly before repeating loop
        plt.pause(0.1)

# Call updatedLine in response to slider updates
slider1.on_changed(updateLine)
slider2.on_changed(updateLine)

# Create start and stop buttons
axButton1 = plt.axes([0.1,0.9,0.1,0.1])
btn1 = Button(axButton1,'Start')
axButton2 = plt.axes([0.25,0.9,0.1,0.1])
btn2 = Button(axButton2,'Stop')

# Define function to set onOff=False and stop continuous plot updating
def stopPlot(val):
    global onOff
    onOff = False

# Set responses to button clicks. Button 1 calls updateLine with current slider values.
btn1.on_clicked(updateLine)
# Button 2 calls function to stop continuous plotting
btn2.on_clicked(stopPlot)

# Show plot in GUI window
plt.show()
