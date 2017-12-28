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
        self.load_experimental_data()
        self.load_theoretical_data()
        self.nice_plots()
        self.plot_raw()
        self.plot_gamma()

    def calc_exp_bins(self):
        channel_numbers = np.arange(1, 1026, 1)
        k = 52
        self.channels = channel_numbers * 10 + k

    def load_experimental_data(self):
        # time correction factor
        c_t = 1 / (10 * 60)  # 1 / (min * s/min)

        # power correction factor
        c_p = 250 / 1E5  # 250 W(th) / 100 kW(th)

        # detector efficiency
        e = 1
        c_e = 1 / e

        # combine into one constant
        c = c_t * c_p * c_e

        spectrum_data = np.loadtxt('nai_gamma_spectrum.tka', skiprows=1) * c
        self.exp_gamma_raw = Spectrum(self.channels, spectrum_data)

        bg_data = np.loadtxt('background.tka', skiprows=1) * c
        self.background = Spectrum(self.channels, bg_data)

        true_data = spectrum_data - bg_data * c
        self.exp_gamma = Spectrum(self.channels, true_data)

    def load_theoretical_data(self):
        # scale56 bin structure
        self.scale56 = np.loadtxt('scale56.txt') * 1E3  # normalize to keV

        # calculate normalization
        tally_area = tally_area = np.pi * (1.27 ** 2)
        tally_area = 1  # for now not considering the tally area
        c = 8.32 / (200 * 1.60218e-13 * tally_area)
        c *= 250  # normalize to 250 W(th)

        # load theoretical gamma data
        data = np.loadtxt('gdata_total.txt')
        data = data.T[1][1:]
        data *= c

        # create spectrum
        self.the_gamma = Spectrum(self.scale56, data)

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
        ax.set_ylabel('Fluence $s^{-1}keV^{-1}$')
        ax.set_yscale('log')
        ax.set_xlim(0, 1E4)
        ax.plot(self.exp_gamma_raw.step_x, self.exp_gamma_raw.step_y, label='Raw')
        ax.plot(self.background.step_x, self.background.step_y, label='Background')
        ax.plot(self.exp_gamma.step_x, self.exp_gamma.step_y, label='Corrected')
        ax.legend()
        fig.savefig('raw.png', dpi=300)

    def plot_gamma(self):
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Energy $keV$')
        ax.set_ylabel('Fluence $s^{-1}keV^{-1}$')
        ax.set_yscale('log')
        ax.set_xlim(0, 1E4)
        ax.set_ylim(1E-5, 4E3)
        ax.plot(self.the_gamma.step_x, self.the_gamma.step_y, label='MCNP')
        ax.plot(self.exp_gamma.step_x, self.exp_gamma.step_y, label='NaI Detector')
        ax.legend()
        fig.savefig('gamma.png', dpi=300)


Plot()
