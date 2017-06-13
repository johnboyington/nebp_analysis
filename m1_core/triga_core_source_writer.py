'''
Writes the SDEF source used in mcnp for the TRIGA reactor core.

Inputs - fission rate values from the KCODE problem for each fuel element segment.
Outputs - SDEF source cards with both fuel element location and axial distribution.

'''
import numpy as np
import matplotlib.pyplot as plt


def cardWriter(card, data, elements):  
    '''
    Function: cardWriter
    
    This will write multiline cards for SI and SP distributions for mcnp inputs
    
    Input Data:
        card - name and number of the card
        data array - a numpy array containing the data you'd like placed in the card.
        Outputs:
            a string that can be copied and pasted into an mcnp input file
    '''
    s = '{}   '.format(card)
    empty_card = '   ' + ' ' * len(card)
    elements_per_row = elements
    row_counter = 0
    element = '{:6}  ' if data.dtype in ['int32', 'int64'] else '{:14.6e}  '
    for i, d in enumerate(data):
        s += element.format(d)
        row_counter += 1
        if row_counter == elements_per_row and i + 1 != len(data):
            row_counter = 0
            s += '\n{}'.format(empty_card)
    s += '\n' 
    return s

#load in rod position data
positions = np.loadtxt('input/pos.txt')

# load in magnitudes from mcnp output
mag_file = open('input/mag.txt', 'r').readlines()
mag = np.array([line.split()[0] for i, line in enumerate(mag_file) if i % 4 == 2]).astype(float).reshape(-1,7)

#calculate the magnitude of each individual rod by summing each segment
rodmag = np.sum(mag, axis= 1)

#try to redo axdiv calculation
axial_divisions = (np.array(range(8)) - 3.5) * (2 * 19.051) / 7.

#create distributed source numbers
position_number = np.array(range(len(positions[:,0]))) + 101


#write the source as a string and write it as a file
s  = 'c ******************************************************************************\n'
s += 'c                         TRIG CORE SOURCE SPECIFICATION                        \n'
s += 'c ******************************************************************************\n'
s += 'SDEF ERG=D1 RAD=D2  AXS=0 0 1  POS=D3  EXT=FPOS=D5 \n'
s += 'SP1  -3\n'
s += 'SI2   0  1.8669\n'
s += 'SP2 -21  1\n'
s += cardWriter('SI3  L', positions.flatten(), 3)
s += cardWriter('SP3   ', rodmag, 4)

s += cardWriter('DS5  S', position_number, 8)


for i, pos in enumerate(position_number):
    s += cardWriter('SI{}  H'.format(pos), axial_divisions, 4)
    s += cardWriter('SP{}  D'.format(pos), np.insert(mag[i], 0, 0), 4)
    
with open('output/core_source.txt', 'w') as H:
            H.write(s)

plt.figure(0, figsize=(4,4))
plt.xlim(-21, 21)
plt.ylim(-21, 21)
plt.scatter(positions[:,0], positions[:,1], s=250)
plt.savefig('output/source_rod_locations.png')