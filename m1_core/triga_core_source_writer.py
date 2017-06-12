#this code will help write the cards for the axial dependant fission rate triga calc

import numpy as np
import matplotlib.pyplot as plt


def cardWriter(card, data):  
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
    elements_per_row = 5
    row_counter = 0
    mystring = '{:6}  ' if data.dtype in ['int32', 'int64'] else '{:.6e}  '
    for d in data:
        s += mystring.format(d)
        row_counter += 1
        if row_counter == elements_per_row:
            row_counter = 0
            s += '\n{}'.format(empty_card)
    s += '\n'
    return s

positions = np.loadtxt('input/pos.txt')
print positions [0]


mag_file = open('input/mag.txt', 'r').readlines()
mag = np.array([line.split()[0] for i, line in enumerate(mag_file) if i % 4 == 2]).astype(float)
print mag

rodmag = []
for ii in range(len(mag)/7):
    summa = 0
    for jj in range(7):
        summa += mag[(7 * ii) + jj]
    rodmag.append(summa)
print len(rodmag)
mag = mag.reshape(-1, 7)

axlen = (2 * 19.051) / 7
axdiv=[]
for kk in range(8):
    axdiv.append(-19.051 + (axlen * kk))


s = ''
s += 'c ******************************************************************************\n'
s += 'c                             SOURCE SPECIFICATION                              \n'
s += 'c ******************************************************************************\n'
s += 'SDEF ERG=D1 RAD=D2  AXS=0 0 1  POS=D3  EXT=FPOS=D5 \n'
s += 'SP1 -3\n'
s += 'SI2   0  1.8669\n'
s += 'SP2 -21  1\n'
s += cardWriter('SI3  L', positions.flatten())
s += cardWriter('SP3   ', np.array(rodmag))

position_number = np.array(range(len(positions[:,0]))) + 100
s += cardWriter('DS5  S', position_number)


for pos in range(len(position_number)):
    s += cardWriter('SI{}  H'.format(pos), np.array(axdiv))
    s += cardWriter('SP{}  D'.format(pos), np.concatenate([0], mag[pos]))
print s
    
with open('output/core_source_altered.txt', 'w') as H:
            H.write(s)

plt.figure(0, figsize=(4,4))
plt.xlim(-21, 21)
plt.ylim(-21, 21)
plt.scatter(xpos, ypos, s=250)
plt.savefig('output/source_rod_locations.png')