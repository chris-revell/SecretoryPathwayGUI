import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button

initA = 1.0
initB = 5.0

# Set up plot canvas
fig,ax=plt.subplots()
plt.subplots_adjust(left=0.1,bottom=0.35)
ax.set_xlim([0,11])
ax.set_ylim([0,10])
# Initialise a flat line
x = np.arange(0,11,0.01)
y = np.zeros(np.shape(x))
line, = ax.plot(x,y)

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

# Function to continually update the plot so long as onOff=True
def updateLine(val):
    global onOff
    # Set onOff = True if currently False
    if onOff:
        pass
    else:
        onOff = True
    # Extract A and B values from sliders
    a_val = slider1.val
    b_val = slider2.val
    while onOff:
        # Shift most of the values left by 20 points
        y[:-20]=y[20:]
        # Fill the last 20 points at the right of the plot with Gaussian noise adjusted with A and B
        y[-20:]=a_val*np.random.randn(20)+b_val
        # Redraw line
        line.set_ydata(y)
        plt.draw()
        # Stop briefly before repeating loop
        plt.pause(1)

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
