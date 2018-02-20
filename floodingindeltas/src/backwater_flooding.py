# This is a module to demonstrate Qw-Qt relation
# The module is written and executed in Python
# This is developed by Kensuke Naito, University of Illinois at Urbana-Champaign
# Email: knaito2@illinois.edu
# Modification for the SedEdu deployment by Andrew J. Moodie (amoodie@rice.edu)


# IMPORT LIBLARIES
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

# SET PARAMETERS
def rootInit():
    L = 1600e3; # length of domain
    nx = 400; # number of nodes
    dx = L / nx; # width of cells
    x = np.arange(0, L+1, dx); # define x-coordinates

    start = 43; # pin-point to start eta from
    S0 = 7e-5; # bed slope
    eta = np.linspace(start, start - S0*(L), nx+1); # channel bed

    mou = 0.75; # fraction of x channelized (i.e. mouth position)
    thet = 2; # plume spreading angle

    B0 = 1100; # basic channel width
    # B = set_B(B0, mou, thet, nx, dx); # channel width
    # print(B)
    # S = get_slope(eta, nx, dx); # bed slope at each node
    H0 = 0; # 
    Cf = 0.0047;
    Qwinit = 10000;
    Qw = Qwinit;
    Qwbf = 35000;
    Qwmax = 60000;
    Qwmin = 5000;

    # [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx);
    # [Xs] = find_backwaterregion(H, dx);
    # [zed] = 0.5 + get_backwater_dBdx(eta, S, B, H0, Cf, Qwbf, nx, dx);


    # MAIN ROUTINE



    # setup the figure
    plt.rcParams['toolbar'] = 'None'
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.15, bottom=0.4)
    l, = plt.plot(x/1000, eta, lw=2, color='blue') # plot initial condition
    # plt.title("")
    ax.set_xlabel("distance from Head of Passes (km)")
    ax.set_ylabel("elevation (m)")
    plt.ylim(-50, 100)
    # plt.grid(b=True, which='both')

    # setup sliders
    ax_color = 'lightgoldenrodyellow'
    ax_Qw = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=ax_color)
    slide_Qw = Slider(ax_Qw, 'water discharge (m^3/s)', Qwmin, Qwmax, valinit=100)

    # setup Reset bottun
    resetax = plt.axes([0.8, 0.01, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=ax_color, hovercolor='0.975')


    # connect sliders
    slide_Qw.on_changed(update)

    # update()
    button.on_clicked(reset)

    # show the results
    plt.show()







def update(val):
    # read values from the sliders


    # conpute sediment load
    print("test")

    # plot results
    # l.set_ydata(Qt)
    # fig.canvas.draw_idle()


def reset(event):
    slide_Qw.reset()


# setup function to update results





# function [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
#     % backwater formulated for changing width
#     H = NaN(1,nx+1); % preallocate depth 
#     H(nx+1) = abs(H0 - eta(nx+1)); % water depth at downstream boundary
#     for i = nx:-1:1
#         % predictor step: computation of a first estimation of the water depth Hp
#         [Frsqp] = get_froude(Qw, H(i+1), B(i)); % calculate froude from conditions at i+1
#         [dBdx] = (B(2:nx+1) - B(1:nx)) ./ dx;
#         [dHdxp] = get_dHdx_dBdx(S(i+1), Cf, Frsqp, H(i+1), B(i), dBdx(i)); % get dHdx width changing form
#         Hp = H(i+1) - dHdxp * dx; % solve for H prediction
#         % corrector step: computation of H
#         [Frsqc] = get_froude(Qw, Hp, B(i)); % calculate froude at i with prediction depth
#         [dHdxc] = get_dHdx_dBdx(S(i), Cf, Frsqc, Hp, B(i), dBdx(i)); % doaa
#         % convolution of prediction and correction, trapezoidal rule
#         H(i) = H(i+1) - ( (0.5) * (dHdxp + dHdxc) * dx );
#     end
    
#     function [dHdx] = get_dHdx_dBdx(S_loc, Cf, Frsq, H_loc, B_loc, dBdx)
#         % formulation to get dHdX for a changing width backwater formulation
#         dHdx = ((S_loc-(Cf*Frsq))/(1-Frsq)) + (Frsq/(1-Frsq)) * (H_loc/B_loc) * (dBdx);
#     end
  
#     function [Frsq] = get_froude(Qw, H, B)
#         g = 9.81; % gravitational acceleration constant
#         Frsq = ( Qw^2 / (g * B^2* H^3) );
#     end
# end

# function [Xs] = find_backwaterregion(H, dx)
#     dHdx = (H(2:end)-H(1:end-1))/dx;
#     dHdxdx = (dHdx(2:end)-dHdx(1:end-1))/dx;
#     threshold = dHdxdx*1e9 < 0.080;
#     changepts = abs(threshold(2:end) - threshold(1:end-1));
#     Xs = find(changepts, 2) * dx;
# end

def set_B(B0, mou, thet, nx, dx):
    B = np.zeros(1, nx+1);
    chanLen = mou * nx;
    B[1:chanLen] = B0;
    B[chanLen+1:end] = B0 + 2*((np.arange(1, (nx+1-chanLen+1))) * dx * np.tan(np.radians(thet)));
    return B

# def get_slope(eta, nx, dx):
#     # return slope of input bed (eta)
#     S = np.zeros(1,nx+1);
#     S[1] = (eta[1] - eta[2]) / dx;
#     S[2:nx] = (eta[1:nx-1] - eta[3:nx+1]) / (2*dx);
#     S(nx + 1) = (eta(nx) - eta(nx + 1)) / dx;
#     return S









if __name__ == '__main__':
    # app = QApplication(sys.argv)
    root = rootInit()
    # root.show()
    # sys.exit(app.exec_())