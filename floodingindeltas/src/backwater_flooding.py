# interactive backwater module written by Andrew J. Moodie
# see module and website for more information
# classroom module for this model can be found at 
# http://www.coastalsustainability.rice.edu/outreach/
# the model setup below is parameterized to the Lower Mississippi River
# as established by Nittrouer et al., 
# Spatial and temporal trends, GSAB, 2012


# IMPORT LIBLARIES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widget
import channel, hydro


# SET PARAMETERS
# def rootInit():
L = 1600e3 # length of domain
nx = 400 # number of nodes
dx = L / nx # width of cells
x = np.arange(0, L+1, dx) # define x-coordinates

start = 43 # pin-point to start eta from
S0 = 7e-5 # bed slope
eta = np.linspace(start, start - S0*(L), nx+1) # channel bed

mou = 0.75 # fraction of x channelized (i.e. mouth position)
thet = 2 # plume spreading angle

B0 = 1100 # basic channel width
B = channel.set_B(B0, mou, thet, nx, dx) # channel width
S = channel.get_slope(eta, nx, dx) # bed slope at each node
H0 = 0 # 
Cf = 0.0047
Qwinit = 10000
Qw = Qwinit
Qwbf = 35000
Qwmax = 60000
Qwmin = 5000

H = hydro.get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
Xs = hydro.find_backwaterregion(H, dx)
zed = 0.5 + hydro.get_backwater_dBdx(eta, S, B, H0, Cf, Qwbf, nx, dx)

# setup the figure
plt.rcParams['toolbar'] = 'None'
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.4, top=0.95, right=0.95)
etaLine, = plt.plot(x/1000, eta, lw=2, color='black') # plot initial condition
waterLine, = plt.plot(x/1000, eta+H, lw=2, color='blue') # plot initial condition
# plt.title("")
ax.set_xlabel("distance from Head of Passes (km)")
ax.set_ylabel("elevation (m)")
plt.ylim(-50, 100)
# plt.grid(b=True, which='both')

# setup sliders
ax_color = 'lightgoldenrodyellow'
ax_Qw = plt.axes([0.25, 0.25, 0.5, 0.05], facecolor=ax_color)
slide_Qw = widget.Slider(ax_Qw, 'water discharge (m^3/s)', Qwmin, Qwmax, valinit=Qwinit, valfmt="%1.0f")

# setup Reset bottun
resetax = plt.axes([0.8, 0.01, 0.1, 0.04])
button = widget.Button(resetax, 'Reset', color=ax_color, hovercolor='0.975')










def update(val):
    # read values from the sliders
    Qw = slide_Qw.val
    H = hydro.get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx);
    
    waterLine.set_ydata(eta+H)
    fig.canvas.draw_idle()



def reset(event):
    slide_Qw.reset()


# connect sliders
slide_Qw.on_changed(update)

# update()
button.on_clicked(reset)

# show the results
plt.show()


# if __name__ == '__main__':
#     # app = QApplication(sys.argv)
#     root = rootInit()
#     # root.show()
#     # sys.exit(app.exec_())