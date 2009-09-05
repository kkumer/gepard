""" 
Auxilliary plotting functions written for ad-hoc plotting
of some specific datasets. Note that actual plotting algorithm
is in `utils` module.

plotHALLA -- plots HALL-A data

"""


import os

import pylab as plt
from matplotlib.ticker import MultipleLocator

import Data
import utils

def plotHERMES(data, approach, pars, path=None, fmt='png'):
    """Makes plot of HERMES preliminary BCA and BSA data with fit line defined by pars"""

    #2: "data/ep2epgamma-BCA-HERMES-08-cos1_b.dat", 
    #4: "data/ep2epgamma-BCA-HERMES-08-cos0_b.dat",
    #5: "data/ep2epgamma-ALU-HERMES-08-Isin1.dat"]
    ids = [2, 4, 5]
    title = 'HERMES-08'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['mt', 'xB', 'Q2']
    ylims = [(-0.05, 0.24), (-0.15, 0.05), (-0.30, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y in range(3):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(MultipleLocator(0.1))  # tickmarks
            utils.subplot(ax, data[ids[y]][x*6:x*6+6], xaxes[x], fits=[(approach, pars)])
            apply(ax.set_ylim, ylims[y])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotCLAS(data, approach, pars, path=None, fmt='png'):
    """Makes plot of CLAS BSA data with fit line defined by pars"""

    #datafile = "data/ep2epgamma-ALU-CLAS_KK-07.dat" # id = 25
    dataset = data[25]
    title = 'CLAS-07'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    nmax = len(dataset) - 1
    # Each different Q2 has its own panel
    panel = 1
    npt = 0
    panelpoints = []
    for panel in [7, 8, 9, 4, 5, 6, 1, 2, 3]:
        ax = fig.add_subplot(3,3,panel)
        panelpoints = []
        pt = dataset[npt]
        Q2 = pt.Q2
        while Q2 == pt.Q2:
            panelpoints.append(pt)
            npt += 1
            if npt == nmax:
                break
            pt = dataset[npt]
        # now transform list into DataSet instance ...
        panelset = Data.DataSet(panelpoints)
        panelset.__dict__ = dataset.__dict__.copy()
        # ... and plot
        utils.subplot(ax, panelset, 'mt', ['Q2', 'xB'], fits=[(approach, pars)])
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotHALLA(data, approach, pars, path=None, fmt='png'):
    """Makes plot of HALL-A data with fit line defined by pars"""

     # "data/ep2epgamma-XLU-HALLA_KK-06-ImCI.dat",
     # "data/ep2epgamma-XUU-HALLA_KK-06-ReCI.dat",
     # "data/ep2epgamma-XUU-HALLA_KK-06-ReCpDCI.dat",
     # "data/ep2epgamma-XLU-HALLA-06-Q2_19_t_028.dat",
     # "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_017.dat",
     # "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_033.dat"
    ids = [26, 27, 28, 15, 17, 24]
    title = 'HALLA-06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    panel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,panel)
        if panel<4:
            utils.subplot(ax, data[id], 't', fits=[(approach, pars)])
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
        else:
            utils.subplot(ax, data[id], 'phi', ['Q2', 't'], fits=[(approach, pars)])
            ax.xaxis.set_major_locator(MultipleLocator(120))
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotfit(pars):
    plotHALLA(pars)
    plotCLAS(pars)
    plotHERMES(pars)
    return

# FIXME: doesn't work!


def plotHALLAphi(pars):
    """Makes plot of HALL-A cross-section data with fit line defined by pars"""
    datafiles = [
             "data/ep2epgamma-XLU-HALLA-06-Q2_15_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_023.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_028.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_033.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_023.dat",
             "data/ep2epgamma-XLU-HALLA-06-Q2_19_t_028.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_033.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_023.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_028.dat",
             "data/ep2epgamma-XLU-HALLA-06-Q2_23_t_033.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_017.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_023.dat",
             #"data/ep2epgamma-XUU-HALLA-06-Q2_23_t_028.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_033.dat"
             ]
    fig = plt.figure()
    panel = 1
    for file in datafiles: 
        datapoints = []
        loadfile(file, datapoints)
        ax = fig.add_subplot(2,3,panel)
        plotpoints(ax, datapoints, pars, 'phi', ['Q2', 't'])
        #plt.ylim(-0.05, 0.24)
        panel += 1
    fig.show()



def minustest(pars):
    """Tabulates A0-r*A1 combination"""
    A0points = []
    A1points = []
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos0_b.dat", A0points)
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos1_b.dat", A1points)
    print "  A0     A1    |    r     |   m_exp  |  m_the"
    print "---------------------------------------------"
    for n in range(len(A0points)):
        p0 = A0points[n]
        p1 = A1points[n]
        print "% 1.3f % 1.3f  |  % 1.3f  |  % 1.3f  | % 1.3f" % (p0.val, p1.val, p1.r, 
                p0.val - p1.r*p1.val, BCA0minusr1(p0, pars))

def pminustest(pars):
    """Tabulates A0-r*A1 combination"""
    A0points = []
    A1points = []
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos0_a.dat", A0points)
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos1_a.dat", A1points)
    print "  A0     A1    |    r     |   m_exp  |  m_the"
    print "---------------------------------------------"
    for n in range(len(A0points)):
        p0 = A0points[n]
        p1 = A1points[n]
        print "% 1.3f % 1.3f  |  % 1.3f  |  % 1.3f  | % 1.3f" % (p0.val, p1.val, p1.r, 
                p0.val - p1.r*p1.val, BCA0minusr1(p0, pars))

