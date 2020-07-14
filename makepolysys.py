#!/usr/bin/env python
#

from pylab import *
import sys

def ktoi(k, kmin):
    return k-kmin+0

def itok(i, kmin):
    return i+kmin-0

def Cpqk(p,q,k):
    if p+q==k: return 1
    else:      return 0

kmin = -2
kmax = 2

allk = arange(kmin, kmax+1)
f    = [0]*len(allk) #zeros(len(allk))

# Set up f vector
f[ktoi(1, kmin)] = 1

polysysvar='polysys'

coeffs=[]
powers=[]

for ki, k in enumerate(allk):
    eqncoeffs = []
    kpowers   = []
    # Add the -k^2 term
    if (k != 0):
        eqncoeffs.append(-k*k)
        upower = [0]*len(allk)
        upower[ktoi(k, kmin)] = 1
        kpowers.append(upower)

    # Add the Cpqk term
    for p in allk:
        q = k-p  # Cpqk
        if ((kmin <= q) and (q <= kmax)):
            upower = [0]*len(allk)
            if (p<q):
                #print('k= %3i p= %3i q= %3i x2'%(k, p, q))
                eqncoeffs.append(2)
                upower[ktoi(p, kmin)] = 1
                upower[ktoi(q, kmin)] = 1
                kpowers.append(upower)
            elif (p==q):
                #print('k= %3i p= %3i q= %3i SQR'%(k, p, q))
                eqncoeffs.append(1)
                upower[ktoi(p, kmin)] = 2
                kpowers.append(upower)

    # Add the f constant term
    if abs(f[ki])>1.0E-16:
        eqncoeffs.append(f[ki])
        upower = [0]*len(allk)
        kpowers.append(upower)

    # Add it to the system
    coeffs.append(eqncoeffs)
    powers.append(kpowers)

print('%% kmin = %02i'%kmin)
print('%% kmax = %02i'%kmax)

print('clear all')
print('addpath("PNLA_MATLAB_OCTAVE");')

# Print out the system
for ic, coeff in enumerate(coeffs):
    sys.stdout.write('%s{%i,1} = [ '%(polysysvar, ic+1))
    for c in coeff:
        sys.stdout.write('%i '%c)
    sys.stdout.write('];\n')

print('')

for ip, power in enumerate(powers):
    sys.stdout.write('%s{%i,2} = [ '%(polysysvar, ip+1))
    for iv, vec in enumerate(power):
        sys.stdout.write(' '.join([str(x) for x in vec]))
        if (iv<len(power)-1 ): sys.stdout.write('; ')
    sys.stdout.write('];\n')

print('[root d c ns check cr digits] = sparf(polysys);')
print('disp(root)')
