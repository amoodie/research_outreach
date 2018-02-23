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
import channel, hydro, utils


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
mouIdx = int(mou*nx)
thet = 2 # plume spreading angle
RKs = np.array([0, 165, 368, 425, 505])
RKidxs = np.int_( (nx*mou) - np.round(RKs*1000/dx) )

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
plt.rcParams['figure.figsize'] = 11, 7
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.075, bottom=0.5, top=0.95, right=0.95)
ax.set_xlabel("distance from Head of Passes (km)")
ax.set_ylabel("elevation (m)")
plt.ylim(-50, 100)
plt.xlim(L/1000*0.25, L/1000-(L/1000*0.125))
# set(hand.ax, 'xTickLabels', cellfun(@num2str, num2cell(abs((cellfun(@str2num, (get(gca, 'XTickLabels')))) - (L/1000*mou))), 'UniformOutput', false))

# add plot elements
etaLine, = plt.plot(x/1000, eta, lw=2, color='black') # plot bed
zedLine = plt.plot(x[:mouIdx]/1000, eta[:mouIdx]+zed[:mouIdx], \
    'k--', lw=1.2) # plot levee
waterLine, = plt.plot(x/1000, eta+H, lw=2, color='blue') # plot initial condition
QwValue = plt.text(0.7, 0.9, "Qw = " + utils.format_number(Qw), transform=ax.transAxes)
BwValue = plt.text(( (Xs[1]-Xs[0])/4 + Xs[0])/1000, 52, \
    "backwater from \n" + "RK " + str(L*mou/1000-round(Xs[0]/1000)) + " to " + str(L*mou/1000-round(Xs[1]/1000)), \
    horizontalalignment="center", backgroundcolor="white")
BwBracket, = plt.plot(np.array([Xs[0], Xs[0], Xs[1], Xs[1]])/1000, np.array([36, 40, 40, 36]), 'k-', lw=1.2)



# add slider
ax_color = 'lightgoldenrodyellow'
ax_Qw = plt.axes([0.075, 0.35, 0.525, 0.05], facecolor=ax_color)
slide_Qw = utils.MinMaxSlider(ax_Qw, 'water discharge (m$^3$/s)', Qwmin, Qwmax, 
    valinit=Qwinit, valstep=500, transform=ax.transAxes)

# add gui table

ax_overTable = plt.axes([0.20, 0.1, 0.5, 0.1], frameon=False, xticks=[], yticks=[])
tabData = [['0', '0', False], ['0', '0', False],
           ['0', '0', False], ['0', '0', False],
           ['0', '0', False]];
tabRowName = ['Head of Passes (RK 0)', 'New Orleans (RK 165)', 'Baton Rouge (RK 368)',
              'St. Francisville (RK 425)', 'Old River Diversion (RK 505)']
tabColName = ['flow depth (m)', 'stage (m)', 'over levee?'];
overTable = plt.table(cellText=tabData, rowLabels=tabRowName,
                      colLabels=tabColName, colWidths=[0.3, 0.3, 0.3],
                      loc="center")
overTable.scale(1, 1.5) # xscale, yscale
[ overTable._cells[(c, 0)]._text.set_text(utils.format_table(HRK)) 
    for c, HRK in zip(np.arange(1,6), H[RKidxs]) ]

# overTable._cells[(2, 1)]._text.set_text("TEST")

# tab.Data(1:end, 2) = (  sprintfc( '%10.1f', (H(RKidxs)) )  ) # insert proper depths
# tab.Data(1:end, 3) = (  sprintfc( '%10.1f', (eta(RKidxs)+H(RKidxs))' )  ); % insert proper stage


# add gui buttons




# setup Reset bottun
resetax = plt.axes([0.8, 0.01, 0.1, 0.04])
button = widget.Button(resetax, 'Reset', color=ax_color, hovercolor='0.975')


def update(val):
    # read values from the sliders
    Qw = slide_Qw.val
    H = hydro.get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
    Xs = hydro.find_backwaterregion(H, dx)
    
    waterLine.set_ydata(eta+H)
    QwValue.set_text("Qw = " + utils.format_number(Qw))
    BwValue.set_text("backwater from \n" + "RK " + str(L*mou/1000-round(Xs[0]/1000)) + \
        " to " + str(L*mou/1000-round(Xs[1]/1000)))
    BwValue.set_x(((Xs[1]-Xs[0])/4 + Xs[0])/1000)
    BwBracket.set_xdata(np.array([Xs[0], Xs[0], Xs[1], Xs[1]])/1000)

    # updateTable()
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