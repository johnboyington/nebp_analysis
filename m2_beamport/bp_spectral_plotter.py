#this code will plot the relative flux of the angular dist. for each energy group

import numpy as np
import matplotlib.pyplot as plt


#makes plt.plot data out of step data
def makeStep(x,y):
    #assert len(x) - 1== len(y)
    Y = np.array([[yy,yy] for yy in np.array(y)]).flatten()
    X = np.array([[xx,xx] for xx in np.array(x)]).flatten()[1:-1]
    return X,Y


#input power level
power_level = 250 #W


#load neutron simulation data
n_data = np.loadtxt('input/ndata.txt')
g_data = np.loadtxt('input/gdata.txt')
n_data_total = np.loadtxt('input/ndata_total.txt')[1:]
g_data_total = np.loadtxt('input/gdata_total.txt')[1:]
n_cos_total = np.loadtxt('input/ncos_total.txt')[1:]
g_cos_total = np.loadtxt('input/gcos_total.txt')[1:]
n_filter_data = np.loadtxt('input/filtered_ndata.txt')[1:]
n_filter_data[:,1] = n_filter_data[:,1] * 7.53942E-8


#determine how many energy groups and cosine groups there are
group_counter = 0
for ii in range(len(n_data[:,0])):
    group_counter += 1
    if n_data[ii,0] > n_data[ii + 1,0]: 
        number_of_erg_groups = group_counter
        group_counter = 0
        break
print 'This data has {} energy groups.'.format(number_of_erg_groups - 1)
number_of_cos_groups = len(n_data[:,0]) / number_of_erg_groups
print 'This data has {} cosine groups.'.format(number_of_cos_groups - 1)


#determine how many energy groups and cosine groups there are (gamma)
group_counter = 0
for ii in range(len(g_data[:,0])):
    group_counter += 1
    if g_data[ii,0] > g_data[ii + 1,0]: 
        g_number_of_erg_groups = group_counter
        group_counter = 0
        break
print 'This data has {} energy groups.'.format(g_number_of_erg_groups - 1)
g_number_of_cos_groups = len(g_data[:,0]) / g_number_of_erg_groups
print 'This data has {} cosine groups.'.format(g_number_of_cos_groups - 1)


#energy groups
erg_groups = n_data[0:number_of_erg_groups,0]
erg_group_width = erg_groups[1:] - erg_groups[:-1]
erg_group_midpoints = erg_groups[:-1] - erg_groups[1:]
erg_group_midpoints = np.append(erg_groups[1:] + (erg_group_midpoints / 2.), [0])
print len(erg_group_midpoints), len(n_data_total)

#energy groups
g_erg_groups = g_data[0:g_number_of_erg_groups,0]
g_erg_group_width = g_erg_groups[1:] - g_erg_groups[:-1]

#neutron flux & error data stored in [angle,group] format
n_flux = n_data[:,1].reshape(number_of_cos_groups,number_of_erg_groups)[1:,1:]
n_error = n_data[:,2].reshape(number_of_cos_groups,number_of_erg_groups)[1:,1:]

#gamma flux & error data stored in [angle,group] format
g_flux = g_data[:,1].reshape(g_number_of_cos_groups, g_number_of_erg_groups)[1:,1:]
g_error = g_data[:,2].reshape(g_number_of_cos_groups, g_number_of_erg_groups)[1:,1:]

#cosine groups
cos_bins = np.loadtxt('input/neutron_cos_bin_structure.txt')
cos_bin_width = cos_bins[:-1] - cos_bins[1:]
cos_midpoints = np.append(cos_bins[1:] + (cos_bin_width / 2.), [0])
assert len(cos_bins) == number_of_cos_groups

#cosine groups (gamma)
g_cos_bins = np.loadtxt('input/gamma_cos_bin_structure.txt')
g_cos_bin_width = g_cos_bins[:-1] - g_cos_bins[1:]
g_cos_midpoints = np.append(g_cos_bins[1:] + (g_cos_bin_width / 2.), [0])
assert len(g_cos_bins) == g_number_of_cos_groups

#create a normalization matrix and normalize the flux values
width_matrix = np.outer(cos_bin_width, erg_group_width)
g_width_matrix = np.outer(g_cos_bin_width, g_erg_group_width)
n_flux_norm = n_flux / width_matrix
g_flux_norm = g_flux / g_width_matrix

#calculate power normalization coefficient and multiply matrices by that
tally_area = np.pi * (1.27 ** 2)
response_thermal = (2.057E5 * 3.63376E-4 ) / (6.242E12) #g * response / (J/Mev)
n_response_tally = 7.53942E-8
g_response_tally = 6.90356E-8

#calc f/t ratio
print 'Sum1 {}'.format(np.sum(n_data_total[0:47,1]))
print 'Sum2 {}'.format(np.sum(n_data_total[47:,1]))
print 'The F/T Ratio is {}'.format(np.sum(n_data_total[0:47,1]) / np.sum(n_data_total[:,1]))


print 'The energy integrated neutron flux is {:7.5e} cm-2s-2 at {} W'.format((n_response_tally / (response_thermal * tally_area)) * power_level, power_level)
print 'The energy integrated gamma flux is {:7.5e} cm-2s-2 at {} W'.format(((g_response_tally * (8.3 / 2.54)) / (response_thermal * tally_area)) * power_level, power_level)

Cn = 1 / (response_thermal * tally_area)
Cg = 1 / (response_thermal * tally_area)
Cg *= (8.3 / 2.54)

print n_data_total[:,1]
print Cn

n_flux_norm *= Cn
g_flux_norm *= Cg
n_data_total[:,1] *= Cn
g_data_total[:,1] *= Cg
n_filter_data[:,1] *= Cn

n_flux_norm *= power_level
g_flux_norm *= power_level
n_data_total[:,1] *= power_level
g_data_total[:,1] *= power_level
n_filter_data[:,1] *= power_level





#reverse the cosine bins for the sake of plotting
cos_bins = cos_bins[::-1]
g_cos_bins = g_cos_bins[::-1]
###############################################################################
#      PLOT 0
#               spectral comparison
###############################################################################
plt.figure(30, figsize=(6, 4.25))
plt.xlabel('Energy  $MeV$')
plt.ylabel('$\phi(E)$  $cm^{-2}s^{-1}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-8, 2E1)
plt.xticks(np.logspace(-8, 1, 10))
plt.tick_params(axis='both', which='minor', bottom='off', top='off', left='off', right='off', labelbottom='off')
plt.plot(makeStep(erg_groups, n_data_total[:,1])[0], 
         makeStep(erg_groups, n_data_total[:,1] / erg_group_width)[1],
         label='Source Spectrum', color='mediumblue')
plt.errorbar(erg_group_midpoints[:-1], 
             n_data_total[:,1] / erg_group_width, 
             (n_data_total[:,1] / erg_group_width) * n_data_total[:,2], 
             color='mediumblue', linestyle="None", capsize=0)
plt.plot(makeStep(erg_groups, n_filter_data[:,1])[0], 
         makeStep(erg_groups, n_filter_data[:,1] / erg_group_width)[1],
         label='Filtered Spectrum', color='mediumseagreen')
plt.errorbar(erg_group_midpoints[:-1], 
             n_filter_data[:,1] / erg_group_width, 
             (n_filter_data[:,1] / erg_group_width) * n_filter_data[:,2], 
             color='mediumseagreen', linestyle="None", capsize=0)
plt.legend()
plt.savefig('output_plot/spectral_comparison.png', dpi=200)
plt.close()

###############################################################################
#      PLOT 0.5
#               spectral comparison (lethargy)
###############################################################################

plt.figure(31, figsize=(6, 4.25))
plt.xlabel('Energy  $MeV$')
plt.ylabel('$\phi(E) E$  $MeV cm^{-2}s^{-1}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-8, 2E1)
plt.xticks(np.logspace(-8, 1, 10))
plt.tick_params(axis='both', which='minor', bottom='off', top='off', left='off', right='off', labelbottom='off')
plt.plot(makeStep(erg_groups, n_data_total[:,1])[0], 
         makeStep(erg_groups, erg_group_midpoints[:-1] * n_data_total[:,1] / erg_group_width)[1],
         label='Source Spectrum', color='mediumblue')
plt.errorbar(erg_group_midpoints[:-1], 
              erg_group_midpoints[:-1] * n_data_total[:,1] / erg_group_width, 
             (n_data_total[:,1] / erg_group_width) * n_data_total[:,2] * erg_group_midpoints[:-1], 
             color='mediumblue', linestyle="None", capsize=0)
plt.plot(makeStep(erg_groups, n_filter_data[:,1])[0], 
         makeStep(erg_groups,  erg_group_midpoints[:-1] * n_filter_data[:,1] / erg_group_width)[1],
         label='Filtered Spectrum', color='mediumseagreen')
plt.errorbar(erg_group_midpoints[:-1], 
              erg_group_midpoints[:-1] * n_filter_data[:,1] / erg_group_width, 
             (n_filter_data[:,1] / erg_group_width) * n_filter_data[:,2] *  erg_group_midpoints[:-1], 
             color='mediumseagreen', linestyle="None", capsize=0)
plt.legend()
plt.savefig('output_plot/spectral_comparison_lethargy.png', dpi=200)
plt.close()





###############################################################################
#      PLOT 1
#               plot flux as a function of energy integrated over angle
###############################################################################
plt.figure(0, figsize=(6, 4.25))
plt.xlabel('Energy  $MeV$')
plt.ylabel('$\phi(E)$  $cm^{-2}s^{-1}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-8, 2E1)
plt.xticks(np.logspace(-8, 1, 10))
plt.tick_params(axis='both', which='minor', bottom='off', top='off', left='off', right='off', labelbottom='off')
plt.plot(makeStep(erg_groups, n_data_total[:,1])[0], 
         makeStep(erg_groups, n_data_total[:,1] / erg_group_width)[1],
         label='Neutron Spectrum', color='mediumblue')
plt.errorbar(erg_group_midpoints[:-1], 
             n_data_total[:,1] / erg_group_width, 
             (n_data_total[:,1] / erg_group_width) * n_data_total[:,2], 
             color='mediumblue', linestyle="None", capsize=0)
plt.plot(makeStep(g_erg_groups, g_data_total[:,1])[0], 
         makeStep(g_erg_groups, g_data_total[:,1] / g_erg_group_width)[1],
         label='Gamma Spectrum', color='goldenrod')
plt.errorbar(erg_group_midpoints[:-1], 
             g_data_total[:,1] / g_erg_group_width, 
             (g_data_total[:,1] / g_erg_group_width) * g_data_total[:,2], 
             color='goldenrod', linestyle="None", capsize=0)
#plt.plot(makeStep(g_erg_groups, 8.3 * np.sum(g_flux, axis=0))[0], 
#         makeStep(g_erg_groups, 8.3 * np.sum(g_flux, axis=0) / g_erg_group_width)[1],
#         label='Gamma Spectrum', color='goldenrod')
plt.legend()
plt.savefig('output_plot/flux_vs_energy.png', dpi=200)
plt.close()


###############################################################################
#      PLOT 2
#               plot flux as a function of angle integrated over energy
###############################################################################
plt.figure(1, figsize=(6, 4.25))
plt.xlabel('Angle ($deg$)')
plt.ylabel('Probability Density')
#plt.yscale('log')
plt.xlim(0, 10)
plt.plot(makeStep(cos_bins, np.sum(n_flux, axis=1))[0], 
         makeStep(cos_bins, (np.sum(n_flux, axis=1) / (cos_bin_width * np.sum(n_flux)))[::-1])[1],
         label='Neutron Spectrum', color='mediumblue')
plt.errorbar(cos_midpoints[:-1], 
             (np.sum(n_flux, axis=1) / (cos_bin_width * np.sum(n_flux))),
             (np.sum(n_flux, axis=1) / (cos_bin_width * np.sum(n_flux))) * n_cos_total[:,1],
             color='mediumblue', linestyle="None", capsize=0)
plt.plot(makeStep(g_cos_bins, np.sum(g_flux, axis=1))[0], 
         makeStep(g_cos_bins, (np.sum(g_flux, axis=1) / (g_cos_bin_width * np.sum(g_flux)))[::-1])[1],
         label='Gamma Spectrum', color='goldenrod')
plt.errorbar(cos_midpoints[:-1], 
             (np.sum(g_flux, axis=1) / (cos_bin_width * np.sum(g_flux))),
             (np.sum(g_flux, axis=1) / (cos_bin_width * np.sum(g_flux))) * g_cos_total[:,1],
             color='goldenrod', linestyle="None", capsize=0)
plt.legend()
plt.savefig('output_plot/flux_vs_angle.png', dpi=200)
plt.close()



###############################################################################
#      PLOT 3
#               plot the relative errors associated with each bin
###############################################################################
plt.figure(2)
plt.xlabel('Energy Group')
plt.ylabel('Cosine Group (0 is Normal)')
plt.title('Relative Neutron Error Associated with Each Bin')
plt.imshow(n_error[::-1,::-1], cmap='plasma', vmin=0, vmax=1.0)
plt.colorbar()
plt.savefig('output_plot/neutron_bin_error.png')
plt.close()


###############################################################################
#      PLOT 4
#               plot the relative errors associated with each bin relative to the maximum error
###############################################################################
plt.figure(3)
plt.xlabel('Energy Group')
plt.ylabel('Cosine Group (0 is Normal)')
plt.title('Relative Neutron Error Associated with Each Bin')
plt.imshow(n_error[::-1,::-1], cmap='plasma', vmin=0, vmax=np.max(n_error))
plt.colorbar()
plt.savefig('output_plot/neutron_bin_error_relative.png')
plt.close()

###############################################################################
#      PLOT 5
#               plot the relative errors associated with each bin (gamma)
###############################################################################
plt.figure(4)
plt.xlabel('Energy Group')
plt.ylabel('Cosine Group (0 is Normal)')
plt.title('Relative Gamma Error Associated with Each Bin')
plt.imshow(g_error[::-1,::-1], cmap='plasma', vmin=0, vmax=1.0)
plt.colorbar()
plt.savefig('output_plot/gamma_bin_error.png')
plt.close()


###############################################################################
#      PLOT 6
#               plot the relative errors associated with each bin relative to the maximum error (gamma)
###############################################################################
plt.figure(5)
plt.xlabel('Energy Group')
plt.ylabel('Cosine Group (0 is Normal)')
plt.title('Relative Gamma Error Associated with Each Bin')
plt.imshow(g_error[::-1,::-1], cmap='plasma', vmin=0, vmax=np.max(g_error))
plt.colorbar()
plt.savefig('output_plot/gamma_bin_error_relative.png')
plt.close()

###############################################################################
#      Error Statistics
#               print the number of bins with each type of error
###############################################################################
error_groups = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 1.0]
error_counter = np.zeros(len(error_groups) - 1)

for e in n_error.flatten():
    for e_group in range(len(error_groups) - 1):
        if e > error_groups[e_group] and e < error_groups[e_group + 1]:
            error_counter[e_group] += 1
error_percentages = (error_counter / np.sum(error_counter)) * 100.

error_string = 'Neutron Bin Error Information\n'
for ii in range(len(error_counter)):
    error_string += '{:3d} bins ({:2.2f}% of total) have a relative error between {} and {}\n'.format(int(error_counter[ii]), error_percentages[ii], error_groups[ii], error_groups[ii + 1])
with open('output_plot/neutron_error_information.txt', 'w') as F:
    F.write(error_string)

#for gammas
g_error_groups = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 1.0]
g_error_counter = np.zeros(len(g_error_groups) - 1)

for e in g_error.flatten():
    for e_group in range(len(g_error_groups) - 1):
        if e > g_error_groups[e_group] and e < g_error_groups[e_group + 1]:
            g_error_counter[e_group] += 1
g_error_percentages = (g_error_counter / np.sum(g_error_counter)) * 100.

g_error_string = 'Gamma Bin Error Information\n'
for ii in range(len(g_error_counter)):
    g_error_string += '{:3d} bins ({:2.2f}% of total) have a relative error between {} and {}\n'.format(int(g_error_counter[ii]), g_error_percentages[ii], g_error_groups[ii], g_error_groups[ii + 1])
with open('output_plot/gamma_error_information.txt', 'w') as F:
    F.write(g_error_string)

###############################################################################
#      Ratios Statistics
#               print the g/n and f/t ratio for each group
###############################################################################

###############################################################################
#      PLOTS 4 and 100-132
#               plot flux as a function of energy integrated over angle
###############################################################################
#for e in range(number_of_erg_groups):
#    plt.figure(100 + e)
#    plt.xlabel('Angle $deg$')
#    plt.ylabel('Relative Flux')
#    plt.title('Angular PDF of Energy Group {}'.format(number_of_erg_groups - (e + 1)))
#    plt.xlim(0, 90)
#    plt.yscale('log')
#    plt.plot(makeStep(cos_bins, n_flux[:,e])[0], 
#             makeStep(cos_bins, (n_flux[:,e] / (cos_bin_width * np.sum(n_flux[:,e])))[::-1])[1])
#    plt.errorbar(cos_midpoints[0:number_of_cos_groups - 1], 
#                 n_flux[:,e] / (cos_bin_width * np.sum(n_flux[:,e])), 
#                 (n_flux[:,e] / np.sum(n_flux[:,e])) * n_error[:,e],
#                 linestyle="None", capsize=0)
#    plt.savefig('output_plot/angular_dist_for_each_energy_group/angle_dist_group_{}.png'.format(number_of_erg_groups - (e + 1)))
#
##trying to make a super subplot
#plt.figure(4, figsize=(12,55))
#plt.xlim(0, 90)
#for e in range(number_of_erg_groups - 1):
#    plt.subplot(number_of_erg_groups - 1, 1, e + 1)
#    plt.xlabel('Angle $deg$')
#    plt.ylabel('Relative Flux')
#    plt.title(' Energy Group {}'.format(number_of_erg_groups - (e + 1)))
#    plt.ylim(1E-4, 1)
#    plt.yscale('log')
#    plt.plot(makeStep(cos_bins, n_flux[:,e])[0], 
#             makeStep(cos_bins, (n_flux[:,e] / (cos_bin_width * np.sum(n_flux[:,e])))[::-1])[1])
#    plt.errorbar(cos_midpoints[0:number_of_cos_groups - 1], 
#                 n_flux[:,e] / (cos_bin_width * np.sum(n_flux[:,e])), 
#                 (n_flux[:,e] / np.sum(n_flux[:,e])) * n_error[:,e],
#                 linestyle="None", capsize=0)
#plt.savefig('output_plot/angle_dist.pdf')
###############################################################################