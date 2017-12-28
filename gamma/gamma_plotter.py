'''
This code compares Aaron Hellinger's results for the gamma spectrum departing
from the NEBP.
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from spectrum import Spectrum


class Plot(object):
    def __init__(self):
        self.calc_exp_bins()
        self.load_data()
        self.nice_plots()
        self.plot_raw()

    def calc_exp_bins(self):
        channel_numbers = np.arange(1, 1026, 1)
        k = 52
        self.channels = channel_numbers * 10 + k

    def load_data(self):
        data = np.loadtxt('nai_gamma_spectrum.tka', skiprows=1)
        self.exp_gamma_raw = Spectrum(self.channels, data)
        data = np.loadtxt('background.tka', skiprows=1)
        self.background = Spectrum(self.channels, data)

    def nice_plots(self):
        rc('font', **{'family': 'serif'})
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        rcParams['xtick.labelsize'] = 12
        rcParams['ytick.labelsize'] = 12
        rcParams['lines.linewidth'] = 1.85
        rcParams['axes.labelsize'] = 15
        rcParams.update({'figure.autolayout': True})

    def plot_raw(self):
        fig = plt.figure(0)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Energy $keV$')
        ax.set_yscale('log')
        ax.plot(self.exp_gamma_raw.step_x, self.exp_gamma_raw.step_y)


Plot()