#this code will help write the cards for the axial dependant fission rate triga calc

import numpy as np
import matplotlib.pyplot as plt

F = open('input/pos.txt', 'r').readlines()
stuff = [i[:-1] for i in F]
things = np.array([(j[16:26], j[27:37]) for j in stuff]).astype(float)
xpos = things[:,0]
ypos = things[:,1]


G = open('input/mag.txt', 'r').readlines()
stuff = [i[:-1] for i in G]
doug = np.array([(j[17:28]) for j in stuff])
ronald = np.array([(k[7:13]) for k in stuff])

mag = []
for nn in range(len(doug)):
    if ((nn+2) % 4 == 0):
        mag.append(doug[nn])
mag = np.array(mag).astype(float) * 10000

cel = []
for mm in range(len(ronald)):
    if ((mm) % 4 == 0):
        np.array(cel.append(ronald[mm])).astype(float)
        
rodmag = []
for ii in range(len(mag)/7):
    summa = 0
    for jj in range(7):
        summa += mag[(7 * ii) + jj]
    rodmag.append(summa)

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

card = ''

for xx in range(len(xpos) / 2):
    if (xx == 0):
        card = 'SI3  L'
    else:
        card = '      '
    s += '{}   {:10.6f} {:10.6f}  0.0000     {:10.6f} {:10.6f}  0.0000\n'.format(card, xpos[xx * 2], ypos[xx*2], xpos[(xx*2) + 1], ypos[(xx*2) + 1])
s += '{}   {:10.6f} {:10.6f}  0.0000\n'.format(card, xpos[len(xpos) - 1], ypos[len(ypos) - 1])

for yy in range(len(xpos) / 5):
    if (yy == 0):
        card = 'SP3   '
    else:
        card = '      '
    s += '{}   {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n'.format(card, rodmag[yy * 5], rodmag[yy * 5 + 1], rodmag[yy * 5 + 2], rodmag[yy * 5 + 3], rodmag[yy * 5 + 4])

num = np.array(range(85)) + 100

for zz in range(len(num) / 5):
    if (zz == 0):
        card = 'DS5  S'
    else:
        card = '      '
    s += '{}   {}  {}  {}  {}  {}\n'.format(card, num[zz * 5], num[zz * 5 + 1], num[zz * 5 + 2], num[zz * 5 + 3], num[zz * 5 + 4])


for aa in range(len(num)):
    card = 'SI{}  H'.format(num[aa])
    s += '{}   {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n'.format(card, axdiv[0], axdiv[1], axdiv[2], axdiv[3])
    s += '           {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n'.format(axdiv[4], axdiv[5], axdiv[6], axdiv[7])
    card = 'SP{}  D'.format(num[aa])
    s += '{}     0.000000 {:10.6f} {:10.6f} {:10.6f}\n'.format(card, mag[aa * 7 + 0], mag[aa * 7 + 1], mag[aa * 7 + 2])
    s += '           {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n'.format(mag[aa * 7 + 3], mag[aa * 7 + 4], mag[aa * 7 + 5], mag[aa * 7 + 6])
    
with open('output/core_source.txt', 'w') as H:
            H.write(s)

plt.figure(0, figsize=(4,4))
plt.xlim(-21, 21)
plt.ylim(-21, 21)
plt.scatter(xpos, ypos, s=250)
plt.savefig('output/source_rod_locations.png')