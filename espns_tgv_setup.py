#!/usr/bin/env python
#

import math
import numpy as np
import sys
import argparse

try:
    from scipy.integrate import solve_ivp
    use_odeint = False
    print("%Loading solve_ivp")
except:
    from scipy.integrate import odeint
    use_odeint = True
    print("%Loading odeint")

# Define NchooseK
def nchoosek(n,k):
    """
    Calculates the binomial coefficients
    / n \      n!
    |   | = -------  for n-k >= 0
    \ k /   k!(n-k)!
    """
    if (n-k>=0):
        return math.factorial(n)/math.factorial(k)/math.factorial(n-k)
    else:
        return 0

def Cpqk_F(p, q, k):
    """
    Returns
              / 1  if p+q  = k
      Cpq,k = |
              \ 0  if p+q != k
    """
    if (p+q == k): return 1
    else:          return 0
    
def Cpqk_L(r,s,t):
    """
    Combination function for Laguerre Polynomials
    """
    sum=0
    for p in range(r+1):
        for q in range(s+1):
            coeff=(-1)**(p+q+t)
            sum = sum + coeff*nchoosek(r,p)*nchoosek(s,q)*nchoosek(p+q,p)*nchoosek(p+q,t)
    return sum

def Cpqk_TOT(pfull, qfull, kfull, skip=0):
    p = pfull[skip:]
    q = qfull[skip:]
    k = kfull[skip:]
    Ndim = len(p)-1
    #print("Ndim = "+repr(Ndim))
    Cpqk = Cpqk_L(p[0], q[0], k[0])
    #if (np.abs(Cpqk)>0): print("Cpqk_L = "+repr(Cpqk)+" p= "+repr(p)+" q= "+repr(q)+" k= "+repr(k))
    for i in range(Ndim):
        Cpqk = Cpqk*Cpqk_F(p[1+i+skip], q[1+i+skip], k[1+i+skip])
    return Cpqk

class kindexing:
    """
    Handles all computations related to indexing
    """
    def __init__(self, Ndim, maxk0, kmin, kmax, removeall0=False):
        self.Ndim = Ndim
        self.removeall0 = removeall0
        kminvec = []
        kmaxvec = []
        Nkvec   = []
        # Construct the max and max vectors
        kminvec.append(0)
        kmaxvec.append(maxk0)
        Nkvec.append(maxk0+1)
        zerobit = 0 if removeall0 else 1
        for i in range(Ndim):
            kminvec.append(kmin)
            kmaxvec.append(kmax)
            Nkvec.append(kmax-kmin+zerobit)   # +0/+1 include zero
        self.kmin = kminvec
        self.kmax = kmaxvec
        self.Nkvec = Nkvec
        kvectemp = list(range(kmin,kmax+1))
        if removeall0:
            try:  # +0/+1 include zero
                kvectemp.remove(0)
            except:
                pass
        self.kvec  = kvectemp
        # kdict maps -kmin, ..,0,.., kmax -> 0, 1, 2.. N
        kdict = []
        for i, ki in enumerate(self.kvec): kdict.append((ki, i))
        self.kvecdict = dict(kdict)
        # Construct the i to k index map
        self.i2kmap = []
        for ivar in range(Ndim):
            for iLag in range(self.kmax[0]+1):
                for k1 in self.kvec:
                    for k2 in self.kvec:
                        if Ndim>2:
                            for k3 in self.kvec:
                                indexvec = [ivar, iLag, k1, k2, k3]
                                if ([k1, k2, k3] == [0,0,0]): continue  # +0/+1 include
                                self.i2kmap.append(indexvec)
                        else:
                            indexvec = [ivar, iLag, k1, k2]
                            if  ([k1,k2] == [0,0]): continue  # +0/+1 include
                            self.i2kmap.append(indexvec)
        self.Nlength  = len(self.i2kmap)
        # Construct list of k vectors (no direction)
        self.kveclist = [ x[1:] for x in self.i2kmap[:(np.prod(np.array(self.Nkvec))-1)] ]
        # Construct the reverse hash-map: k(hash) to i
        khashlist=[]
        for i, k in enumerate(self.i2kmap):
            khashlist.append((self.computekhash(k), i))
            #print(' khash %i -> i %i'%(khash,i))
        self.khashdict = dict(khashlist)

    def computekhash(self, k):
        Nklength = []
        for i in range(len(self.Nkvec)+1):
            Nklength.append(np.prod(np.array(self.Nkvec[i:])))
        iresult  = k[0]*Nklength[0] + k[1]*Nklength[1]
        for idim in range(2, 2+self.Ndim-1):
            ik   = self.kvecdict[k[idim]]
            iresult = iresult + ik*Nklength[idim]
        iresult = iresult + self.kvecdict[k[-1]]
        return iresult
        

    def i2k(self, i):
        return self.i2kmap[i]

    def k2i(self, k):
        """ 
        Map a k vector to index i
        """
        khash = self.computekhash(k)
        return self.khashdict[khash]+0

    def checkmapping(self, verbose=False):
        totaldiffi = 0
        for i, k in enumerate(self.i2kmap):
            diffi = i-self.k2i(self.i2k(i))
            totaldiffi += np.abs(diffi)
            if verbose: print(repr(i)+' %-15s '%repr(self.i2k(i))+'\t'+repr(diffi))
        return totaldiffi
    
    def dumpinfo(self):
        print('kmin     = '+repr(self.kmin))
        print('kmax     = '+repr(self.kmax))
        print('kvec     = '+repr(self.kvec))
        print('kvecdict = '+repr(self.kvecdict))
        print('Nkvec    = '+repr(self.Nkvec))
        print('kveclist = '+repr(self.kveclist))
        print('Nlength  = '+repr(self.Nlength))
        return

def ExactTGV2DSolution(kind, k, nu, timedependence=True, deriv=False):
    """
    Returns the exact Taylor-Green Vortex solution in 2D:
    u =  sin(x)*cos(y)*exp(-2*nu*t)
    v = -cos(x)*sin*y)*exp(-2*nu*t)
    """
    if len(k)>4:
        raise NameError("ExactTGV2DSolution only works in 2D")
    idim = k[0]   # Direction, should be 0 or 1
    k0   = k[1]   # Laguerre frequency
    k1   = k[2]   # kx
    k2   = k[3]   # ky
    # zero coefficients for k1 or k2 outside of -1,+1
    if ((abs(k1)>1) or (abs(k2)>1)):
        return 0
        
        ik1 = np.sign(k1)+1  # +0/+1 include zero
    ik2 = np.sign(k2)+1  # +0/+1 include zero
    #print("--> ik1 = "+repr(ik1)+" ik2 = "+repr(ik2))
    km  = [np.sign(k1), np.sign(k2)]
    
    expnut = ((2.0*nu)**k0)/((2.0*nu+1.0)**(k0+1))
    if deriv:  # Set to the derivative
        expnut = (2**k0)*(k0-2.0*nu)*(nu**(k0-1.0))/(2.0*nu+1.0)**(k0+2.0)
#        expnut = ((2.0*k0*(2.0*nu)**(k0-1))/(1.0+2.0*nu)**(k0+1))
#        - ((2.0*(k0+1.0)*(2.0*nu)**k0)/(1.0+2.0*nu)**(k0+2))

    iI = 1j

    if timedependence: Ft = expnut
    else:              Ft = 1.0

    uhat = ((-1)**(idim+1))*0.25*iI*float(km[idim])*Ft
    return uhat

def ExactTGV2DSolutionAllk(kind, nu, timedependence=True, deriv=False):
    if kind.Ndim!=2:
        raise NameError("ExactTGV2DSolutionAllk only works in 2D")
    allUhat = []
    for ks in kind.i2kmap:
        allUhat.append(ExactTGV2DSolution(kind, ks, nu, timedependence=timedependence, deriv=deriv))
    return allUhat

def DiffLaguerreVec(kind, ks):
    Dvec = kind.Nlength*[0]
    for k0step in range(ks[1]+1):
        k0new = ks[1]-k0step
        kvecnew = ks[:]
        kvecnew[1] = k0new
        inew = kind.k2i(kvecnew)
        Dvec[inew] = 1.0
    return Dvec

def MakeLinear(kind, ks, nu):
    i = kind.k2i(ks)
    kappa2 = sum([x*x for x in ks[-kind.Ndim:]])
    #kappa=np.linalg.norm(np.array(k[-kind.Ndim:]))
    L = np.array(DiffLaguerreVec(kind, ks))
    #print(" kappa2 = "+repr(kappa2)+" nu = "+repr(nu))
    L[i] += float(kappa2)*nu
    return np.array(L)

def MakeAdvectionMat(kind, k):
    N    = kind.Nlength
    kdir = k[0]
    Mmat = np.zeros((N,N))
    for p in kind.kveclist:
        for q in kind.kveclist:
            upind = kind.k2i([0]+p)
            vpind = kind.k2i([1]+p)
            vel_qind = kind.k2i([kdir]+q)
            q1    = q[1]
            q2    = q[2]
            Cpqk = Cpqk_TOT(p, q, k[1:], skip=0)
            Mmat[upind][vel_qind] += q1*Cpqk
            Mmat[vpind][vel_qind] += q2*Cpqk
    return 1j*Mmat

def MakePressureMat(kind, k):
    N     = kind.Nlength
    Mmat  = np.zeros((N,N))    
    Ndim  = len(k)-2
    kdir  = k[0]
    km    = k[2+kdir]
    kappa = np.linalg.norm(k[2:])
    for p in kind.kveclist:
        for q in kind.kveclist:
            Cpqk = Cpqk_TOT(p, q, k[1:], skip=0)
            if (np.abs(Cpqk)>0):
                for i in range(Ndim):
                    for j in range(Ndim):
                        pi = p[1+i] 
                        qj = q[1+j] 
                        ujp_ind = kind.k2i([j]+p)
                        uiq_ind = kind.k2i([i]+q)
                        Mmat[ujp_ind][uiq_ind] += pi*qj*Cpqk
    return -1j*(km/kappa**2)*Mmat


def MakeZeqn(kind, k, nu, uhat, u0hat, verbose=False):
    """
    Make the Z-equation for the polnomial N-S equation at wavenumber k
    """
    L           = MakeLinear(kind, k, nu)
    LinearTerms = np.dot(L, uhat)-u0hat[kind.k2i(k)]
    Nmat        = MakeAdvectionMat(kind, k) + MakePressureMat(kind,k)
    u           = np.array(uhat)
    NL          = np.dot(u, np.dot(Nmat, u))
    if verbose:
        #print(" Z k = "+repr(k))
        print(" L   = "+repr(LinearTerms))
        print(" NL  = "+repr(NL))
    return LinearTerms + NL        

def checkZeqns(kind, nu, uhat=None, u0hat=None, verbose=False):
    if uhat  is None: uhat  = ExactTGV2DSolutionAllk(kind, nu)
    if u0hat is None: u0hat = ExactTGV2DSolutionAllk(kind, nu,
                                                     timedependence=False)
    sumZ = 0
    for kk in kind.i2kmap:    
        Z = MakeZeqn(kind, kk, nu, uhat, u0hat, verbose=verbose)
        sumZ += np.abs(Z)
        if verbose: print('k= %-15s'%repr(kk)+'  Total = '+repr(Z))
    return sumZ

def MakeHJac(kind, k, nu, uhat, verbose=False):
    JVec   = np.zeros(kind.Nlength)*1j
    L      = MakeLinear(kind, k, nu)
    JVec  += L
    Nmat   = MakeAdvectionMat(kind, k) + MakePressureMat(kind,k)
    #print('JVec shape = '+repr(np.shape(JVec)))
    JVec  += np.dot(Nmat + np.transpose(Nmat), uhat)
    return JVec

def MakeHJacMat(kind, nu, uhat, verbose=False):
    """
    Make the Jacobian matrix d(H_k)/d(u_p)
    """
    #print("nu = %e"%nu)
    HJac   = np.zeros((kind.Nlength, kind.Nlength))*1j
    for ki, k in enumerate(kind.i2kmap):
        JacVec = MakeHJac(kind, k, nu, uhat, verbose=verbose)
        HJac[ki, :] = JacVec
    return HJac

def MakeHnu(kind, uhat, verbose=False):
    """
    Make the vector d(H_k)/d(nu)
    """
    Hnu    = np.zeros(kind.Nlength)*1j
    for ki, k in enumerate(kind.i2kmap):
        idim = k[0]   # Direction, should be 0 or 1
        k0   = k[1]   # Laguerre frequency
        k1   = k[2]   # kx
        k2   = k[3]   # ky
        Hnu[ki] = (k1*k1 + k2*k2)*uhat[ki]
    return Hnu

def DavidenkoDeriv(uhat, t, kind, v):
    #print("v+t = "+repr(v+t))
    HJac     = MakeHJacMat(kind, v+t, uhat)
    #print(np.real(HJac))
    Hnu      = MakeHnu(kind, uhat)
    duhatdnu = -np.dot(np.linalg.inv(HJac), Hnu) 
    #duhatdnu = ExactTGV2DSolutionAllk(kind, v+t, deriv=True)
    return duhatdnu

def checkDavidenkoDeriv(kind, nu, verbose=False):
    uhat  = ExactTGV2DSolutionAllk(kind, nu)
    u0hat = ExactTGV2DSolutionAllk(kind, nu, timedependence=False)
    Duhat = ExactTGV2DSolutionAllk(kind, nu, deriv=True)

    duhatdnu = DavidenkoDeriv(uhat, 0, kind, nu)
    # Test the DavidenkoDeriv function
    totaldiff = 0
    for ki, k in enumerate(kind.i2kmap):
        diff = duhatdnu[ki] - Duhat[ki]
        totaldiff += np.abs(diff)
        if verbose: print('k= %-15s'%repr(k)+' '+repr(duhatdnu[ki])+' '+repr(Duhat[ki])+' diff = '+repr(np.abs(diff)))
    return totaldiff

def odeintz(func, z0, t, **kwargs):
    """An odeint-like function for complex valued differential equations."""
    # See https://stackoverflow.com/questions/19910189/scipy-odeint-with-complex-initial-values

    # Disallow Jacobian-related arguments.
    _unsupported_odeint_args = ['Dfun', 'col_deriv', 'ml', 'mu']
    bad_args = [arg for arg in kwargs if arg in _unsupported_odeint_args]
    if len(bad_args) > 0:
        raise ValueError("The odeint argument %r is not supported by "
                         "odeintz." % (bad_args[0],))

    # Make sure z0 is a numpy array of type np.complex128.
    z0 = np.array(z0, dtype=np.complex128, ndmin=1)

    def realfunc(x, t, *args):
        z = x.view(np.complex128)
        dzdt = func(z, t, *args)
        # func might return a python list, so convert its return
        # value to an array with type np.complex128, and then return
        # a np.float64 view of that array.
        return np.asarray(dzdt, dtype=np.complex128).view(np.float64)

    result = odeint(realfunc, z0.view(np.float64), t, **kwargs)

    if kwargs.get('full_output', False):
        z = result[0].view(np.complex128)
        infodict = result[1]
        return z, infodict
    else:
        z = result.view(np.complex128)
        return z

def integrate(kind, nu0, nu1, uhat, use_odeint):
    """
    Integrate the homotopy
    """
    tvec          = np.linspace(0, nu1-nu0, 41)
    if use_odeint:
        sol, infodict = odeintz(DavidenkoDeriv, uhat, tvec, args=(kind, nu0),
                                full_output=True, atol=1.0e-10,hmax=5e-3)
        nuvec = nu0 + tvec
    else:
        fun  = lambda t, y: DavidenkoDeriv(y,t,kind,nu0)
        soly = solve_ivp(fun, [0, nu1-nu0], uhat, method='RK45',t_eval=tvec)
        nuvec = nu0 + soly.t
        sol  = soly.y.transpose()
    return nuvec, sol

okcheck = lambda sum, tol: 'OK' if sum<tol else 'NOT OK'

def runTests(kind, nu, verbose=False):
    """
    Run some unit tests the functions above
    """
    # Check the indexing functions
    totaldiffi = kind.checkmapping(verbose=verbose)
    print('%CHECK K-INDEXING  = '+repr(totaldiffi)+' [%s]'%okcheck(totaldiffi,1.0E-10))

    # Check that the Z-equations are satisfied
    sumZ=checkZeqns(kind, nu, verbose=verbose)
    print("%%CHECK Z-EQNS      = %e"%sumZ+' [%s]'%okcheck(sumZ, 1.0E-10))

        # Test the DavidenkoDeriv function
    totaldiff = checkDavidenkoDeriv(kind, nu, verbose=verbose)
    print("%%CHECK DERIVE      = %e"%totaldiff+' [%s]'%okcheck(totaldiffi,1.0E-10))

    return

def writeZeqn(kind, k, nu, u0hat, verbose=False):
    """
    Get a single Zk equation
    """
    eqncoeffs   = []
    eqnpowers   = []

    # Process constant term
    Const       = -u0hat[kind.k2i(k)]
    eqncoeffs.append(Const)
    eqnpowers.append([0]*kind.Nlength)

    # Process the linear term
    L           = MakeLinear(kind, k, nu)
    for li, lcoeff in enumerate(L):
        if (lcoeff != 0):
            powers = [0]*kind.Nlength
            powers[li] = 1
            eqncoeffs.append(lcoeff)
            eqnpowers.append(powers)
        
    # Process the non-linear term
    Nmat        = MakeAdvectionMat(kind, k) + MakePressureMat(kind,k)
    if verbose:
        print(np.linalg.norm(MakeAdvectionMat(kind, k)),
              np.linalg.norm(MakePressureMat(kind, k)),
              np.linalg.norm(Nmat))
        print("")
    for p in range(kind.Nlength):
        for q in range(p, kind.Nlength):
            if (Nmat[p,q] + Nmat[q,p] != 0.0):
                powers = [0]*kind.Nlength
                if (p == q):
                    eqncoeffs.append(Nmat[p,q])
                    powers[p] = 2
                else:
                    eqncoeffs.append(Nmat[p,q]+Nmat[q,p])
                    powers[p] = 1
                    powers[q] = 1
                eqnpowers.append(powers)
    
    return eqncoeffs, eqnpowers

def writePolysys(stream=sys.stdout):
    # User defined variables here
    Ndim    = 2      # Number of dimensions
    kmaxlag = 2      # Maximum Laguerre expansion
    kmax    = 1      # Maximum Fourier  expansion
    kmin    = -kmax
    nu0     = 1.0
    verbose    = False
    removeall0 = True #False
    polysysvar = 'polysys'

    header="""
clear all
addpath("PNLA_MATLAB_OCTAVE");
addpath("SuiteSparse/SPQR/MATLAB");
stol=2.5E-15;
"""    
    kindex  = kindexing(Ndim, kmaxlag, kmin, kmax, removeall0=removeall0)

    runTests(kindex, nu0, verbose=False)

    uexact = ExactTGV2DSolutionAllk(kindex, nu0)
    u0hat  = ExactTGV2DSolutionAllk(kindex, nu0, timedependence=False)

    allcoeffs = []
    allpowers = []
    for kk in kindex.i2kmap:    
        eqnc, eqnp = writeZeqn(kindex, kk, nu0, u0hat, verbose=verbose)
        allcoeffs.append(eqnc)
        allpowers.append(eqnp)
        if verbose:
            print("k = "+repr(kk))
            print(eqnc)
            print(eqnp)
            print("")

    stream.write(header)
    stream.write('nu = %f;\n'%nu0)
    # Write the exact solution
    stream.write('uexact = [ ')
    for x in uexact:
        stream.write('%12.10f%+20.18fj, '%(np.real(x), np.imag(x)))
    stream.write('];\n')
    
    for ic, coeff in enumerate(allcoeffs):
        stream.write('%s{%i,1} = [ '%(polysysvar, ic+1))
        for c in coeff:
            if (np.iscomplex(c)):
                stream.write('%f%+fj '%(np.real(c), np.imag(c)))
            else:
                stream.write('%f '%c)
        stream.write('];\n')

    stream.write('\n')

    # Write the index
    stream.write('% kindex map\n');
    stream.write('kindex = [');
    for kk in kindex.i2kmap:    
        stream.write('%i %i %i %i;\n'%(kk[0], kk[1], kk[2], kk[3]))
    stream.write('];\n')
    
    for ip, power in enumerate(allpowers):
        stream.write('%s{%i,2} = [ '%(polysysvar, ip+1))
        for iv, vec in enumerate(power):
            stream.write(' '.join([str(x) for x in vec]))
            if (iv<len(power)-1 ): stream.write('; ')
        stream.write('];\n')

    stream.write('\n')
    stream.write('[root, d, c, ns, check, cr, digits] = qdsparf2(%s, stol);'%polysysvar)
    stream.write('\n')
    stream.write("disp(sprintf('diff = %e',norm(uexact-root)))")
    stream.write('\n')

    dispsol="""
%Print the solution in pretty form
screenstring='%2i % 2i % 2i % 2i % 9.7fi %i*i/%i\\n';
latexstring='$%i$ & $%i$ &  $%i$ & $%i$ & $\\\\texttt{% .7f}\\\\ii$ & ${%i\\\\ii}/{%i}$~\\\\tabularnewline\\n';
for i=1:length(kindex),
    idim = kindex(i,1);
    k0   = kindex(i,2);
    k1   = kindex(i,3);
    k2   = kindex(i,4);
    fprintf(latexstring,...
            idim+1, kindex(i,2), k1, k2,...
            imag(full(root(i))), (-1)^(idim+1)*kindex(i,3+idim)*(2.0*nu)^k0, 4*(1+2*nu)^(k0+1));
end\n"""
    stream.write(dispsol)
    
    return

def demo():
    """
    Demonstrate the problem
    """
    # User defined variables here
    Ndim    = 2      # Number of dimensions
    kmaxlag = 0      # Maximum Laguerre expansion
    kmax    = 2      # Maximum Fourier  expansion
    kmin    = -kmax
    nu0     = 1.0
    nu1     = 3.0
    verbose = False
    removeall0=True #False

    kindex  = kindexing(Ndim, kmaxlag, kmin, kmax, removeall0=removeall0)
    
    if verbose: kindex.dumpinfo()

    runTests(kindex, nu0, verbose=False)
    
    exactuhat  = ExactTGV2DSolutionAllk(kindex, nu0)
    exactuhat1 = ExactTGV2DSolutionAllk(kindex, nu1)
    nuvec, sol = integrate(kindex, nu0, nu1, exactuhat, use_odeint)
    diff = sum(np.abs(sol[-1,:]-exactuhat1))
    print("CHECK (INTEGRATE) = %e"%diff+' [%s]'%okcheck(diff, 1.0E-8))

    # Generate a figure and compare against exact
    from matplotlib import pyplot as plt
    plt.rcParams.update({'font.size': 16})
    kplotvec = [[0,0,-1,-1],[1,1,1,1]] if kmaxlag > 0 else [[0,0,-1,-1]]
    nondim = False

    plt.figure(figsize=(12,5))
    for ik, k in enumerate(kplotvec):
        ploti = kindex.k2i(k)
        startval = ExactTGV2DSolution(kindex, k,nu0) if nondim else 1.0
        ExactNuSol = [ExactTGV2DSolution(kindex, kindex.i2k(ploti),t)/startval for t in nuvec]
        crosslabel = 'Homotopy' if ik==0 else ''
        plt.plot(nuvec,np.abs(sol[:,ploti]/startval), 'k+',
                 markersize=12, markeredgewidth=1.25, label=crosslabel)
        plt.plot(nuvec,np.abs(ExactNuSol), '-', linewidth=2,
                 label='Exact m='+repr(k[0])+' k='+repr(k[1:]))
    plt.legend()
    plt.xlabel(r'Viscosity $\nu$')
    if nondim:    plt.ylabel(r'$|\hat{u}_m^k(\nu)/\hat{u}_m^k(\nu_0)|$')
    else:         plt.ylabel(r'$|\hat{u}_m^k(\nu)|$')
    plt.tight_layout()
    #plt.savefig('NS_TGV_homotopy.png')
    plt.show()

if __name__ == "__main__":
    #demo()

     # Parse arguments
    parser = argparse.ArgumentParser(description="Set up TGV problem")
    parser.add_argument(
        "-o",
        "--outfile",
        help="save to file",
        default='',
        type=str,
    )
    args = parser.parse_args()
    outstream=sys.stdout
    if (len(args.outfile)>0):
        print("% Saving to "+args.outfile)
        outstream=open(args.outfile, 'w')
    writePolysys(stream=outstream)
    #if (len(args.outfile)>0):
    outstream.close()
    
