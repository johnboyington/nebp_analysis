'''
This code writes an SDEF card for a source given output data from the North East Beam Port
It is given input data that contains flux information grouped into energy and cos bins.
'''

import numpy as np

'''
Function: cardWriter

This will write multiline cards for SI and SP distributions for mcnp inputs

Input Data:
    card - name and number of the card
    data array - a numpy array containing the data you'd like placed in the card.
Outputs:
    a string that can be copied and pasted into an mcnp input file
'''
def cardWriter(card, data):    
    s = ''
    empty_card = '   '
    elements_per_row = 5
    row_counter = 0
    for space in range(len(card)):
        empty_card += ' '
    s += '{}   '.format(card)
    for d in range(len(data)):
        if data.dtype == 'int32' or data.dtype == 'int64': s += '{:6}  '.format(data[d])
        else: s += '{:.6e}  '.format(data[d])
        row_counter += 1
        if row_counter == elements_per_row:
            row_counter = 0
            s += '\n{}'.format(empty_card)
    s += '\n'
    return s


#load neutron simulation data
n_data = np.loadtxt('input/ndata.txt')


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

#energy groups
erg_groups = n_data[0:number_of_erg_groups,0]
erg_group_width = erg_groups[1:] - erg_groups[:-1]

#flux & error data stored in [angle,group] format
n_flux = n_data[:,1].reshape(number_of_cos_groups,number_of_erg_groups)[1:,1:]
n_error = n_data[:,2].reshape(number_of_cos_groups,number_of_erg_groups)[1:,1:]

#cosine groups
cos_bins = np.loadtxt('input/cos_bin_structure.txt')
cos_bin_width = cos_bins[:-1] - cos_bins[1:]
cos_midpoints = np.append(cos_bins[1:] + (cos_bin_width / 2.), [0])
assert len(cos_bins) == number_of_cos_groups

#create a normalization matrix and normalize the flux values
width_matrix = np.outer(cos_bin_width, erg_group_width)
n_flux_norm = n_flux / width_matrix

cos_bins_radians = np.cos(cos_bins * (np.pi / 180))



###############################################################################
#                          SOURCE WRITER
###############################################################################
dist_source_numbers = 100 + np.array(range(number_of_erg_groups))
dummy_distribution = np.append(np.zeros(number_of_cos_groups - 1), 1)

s = ''
s += 'c  ---------------------------------------------------------\n'
s += 'c                    SOURCE SPECIFICATIONS\n'
s += 'c  ---------------------------------------------------------\n'
s += 'SDEF POS=0 0 0 AXS=1 0 0 EXT=0 VEC=1 0 0 ERG=D3 DIR=FERG=D5 RAD=D6 PAR=1\n'
s += 'SI6   0  1.27\n'
s += 'SP6 -21  1\n'
s += cardWriter('SI3 H', erg_groups)
s += cardWriter('SP3 D', np.append(0, np.sum(n_flux, axis=0)))
s += cardWriter('DS5 S', dist_source_numbers)
s += cardWriter('SI100 H', cos_bins_radians)
s += cardWriter('SP100 D', dummy_distribution)
for e in range(number_of_erg_groups - 1):
    s += cardWriter('SI{} H'.format(dist_source_numbers[e] + 1), cos_bins_radians)
    s += cardWriter('SP{} D'.format(dist_source_numbers[e] + 1), np.append(0, n_flux[:,e]))

with open('output_sdef/nsource.txt', 'w') as H:
    H.write(s)
###############################################################################