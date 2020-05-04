import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button

initA = 5.0
initB = 1.0

fig,ax=plt.subplots()
plt.subplots_adjust(left=0.1,bottom=0.35)
x = np.arange(0,11,0.01)
y = x
line, = ax.plot(x,y)
ax.set_xlim([0,11])
ax.set_ylim([0,50])

axSlider1 = plt.axes([0.1,0.2,0.8,0.05])
axSlider1.set_xticks([0.0,100.0])
slider1 = Slider(axSlider1,'a',valinit=5.0,valmin=0,valmax=10)
axSlider1.xaxis.set_visible(True)
axSlider1.set_xticks([0,10])

axSlider2 = plt.axes([0.1,0.1,0.8,0.05])
axSlider2.set_xticks([0.0,100.0])
slider2 = Slider(axSlider2,'b',valinit=1.0,valmin=0,valmax=5,color='orange')
axSlider2.xaxis.set_visible(True)
axSlider2.set_xticks([0,5])

def updateLine(val):
    a_val = slider1.val
    b_val = slider2.val
    #x = np.arange(0,11,0.01)
    #ax.cla()
    #ax.plot(x,a_val*x+b_val)
    line.set_ydata(a_val*x+b_val)
    plt.draw()


slider1.on_changed(updateLine)
slider2.on_changed(updateLine)


axButton1 = plt.axes([0.1,0.9,0.1,0.1])
btn1 = Button(axButton1,'Reset')

axButton2 = plt.axes([0.25,0.9,0.2,0.1])
btn2 = Button(axButton2,'Set b=0')

def resetSliders(event):
    slider1.reset()
    slider2.reset()
btn1.on_clicked(resetSliders)

def setValue(val):
    slider2.set_val(0)
btn2.on_clicked(setValue)


plt.show()
