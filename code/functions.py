from __future__ import division
import numpy as np
from fatiando.constants import G, SI2MGAL
from scipy.sparse import diags
from scipy.linalg import toeplitz
from numpy.linalg import inv

def fast_eq(x,y,z,h,shape,data,itmax):
    '''
    Calculates the estimate physical property distribution of
    an equivalent layer that repreduces a gravity disturbance
    data through an iterative method [1].

    [1] SIQUEIRA, F. C., OLIVEIRA JR, V. C., BARBOSA, V. C., 2017,
    Fast iterative equivalent-layer technique for gravity data
    processing: A method grounded on excess mass constraint",
    Geophysics, v. 82, n. 4, pp. G57-G69.

    input
    x, y: numpy array - the x, y coordinates
    of the grid and equivalent layer points.
    z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.
    shape: tuple - grid size.
    data: numpy array - the gravity disturbance.
    potential field at the grid points.

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    gzp: numpy array - the predicted data.
    '''
    assert x.size == y.size == z.size == h.size == data.size, 'x, y,\
    z, h and data must have the same number of elements'
    #assert h.all() > z.all(), 'The equivalent layer must be beneath\
    #the observation points'

    #Diagonal matrix calculation
    N,diagonal_A = diagonal(x,y,shape)

    #Initial estimative
    rho0 = i_rho(data,diagonal_A)
    
    #Complete sensibility matrix
    A = sensibility_matrix(x,y,z,h,N)

    #Fast Equivalent layer loop
    m_new = fast_loop(data,A,diagonal_A,rho0,itmax)

    #Final predicted data
    gzp = A.dot(m_new)
    
    A = np.zeros((N, N), dtype=np.float)

    return m_new, gzp

def fast_eq_toep(x,y,z,h,shape,data,itmax):
    '''
    Calculates the estimate physical property distribution of
    an equivalent layer that repreduces a gravity disturbance
    data through an iterative method [1]. 
    This implementation uses a fast way to calculate the forward
    problem at each iteration taking advantage of the BTTB
    (Block-Toeplitz Toreplitz-Block) structures.

    [1] SIQUEIRA, F. C., OLIVEIRA JR, V. C., BARBOSA, V. C., 2017,
    Fast iterative equivalent-layer technique for gravity data
    processing: A method grounded on excess mass constraint",
    Geophysics, v. 82, n. 4, pp. G57-G69.

    input
    x, y: numpy array - the x, y coordinates
    of the grid and equivalent layer points.
    z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.
    shape: tuple - grid size.
    data: numpy array - the gravity disturbance.
    potential field at the grid points.

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    gzp: numpy array - the predicted data.
    '''
    assert x.size == y.size == z.size == h.size == data.size, 'x, y,\
    z, h and data must have the same number of elements'
    #assert h.all() > z.all(), 'The equivalent layer must be beneath\
	#the observation points'
    
    #Diagonal matrix calculation
    N,diagonal_A = diagonal(x,y,shape)

    #Initial estimative
    rho0 = i_rho(data,diagonal_A)
    
    #Create first line of sensibility matrix
    BTTB = bttb(x,y,z,h)

    #Fast Equivalent layer loop
    m_new = fast_loop_toep(BTTB,shape,N,data,diagonal_A,rho0,itmax)

    #Final predicted data
    gzp = fast_forward_toep(shape,N,m_new,BTTB)

    return m_new, gzp
    
def fast_eq_bccb(x,y,z,h,shape,data,itmax):
    '''
    Calculates the estimate physical property distribution of
    an equivalent layer that repreduces a gravity disturbance
    data through an iterative method [1]. 
    This implementation uses a fast way to calculate the forward
    problem at each iteration taking advantage of the BTTB
    (Block-Toeplitz Toreplitz-Block) structures.

    [1] SIQUEIRA, F. C., OLIVEIRA JR, V. C., BARBOSA, V. C., 2017,
    Fast iterative equivalent-layer technique for gravity data
    processing: A method grounded on excess mass constraint",
    Geophysics, v. 82, n. 4, pp. G57-G69.

    input
    x, y: numpy array - the x, y coordinates
    of the grid and equivalent layer points.
    z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.
    shape: tuple - grid size.
    data: numpy array - the gravity disturbance.
    potential field at the grid points.

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    gzp: numpy array - the predicted data.
    '''
    assert x.size == y.size == z.size == h.size == data.size, 'x, y,\
    z, h and data must have the same number of elements'
    #assert h.all() > z.all(), 'The equivalent layer must be beneath\
	#the observation points'
    
    #Diagonal matrix calculation
    N,diagonal_A = diagonal(x,y,shape)

    #Initial estimative
    rho0 = i_rho(data,diagonal_A)
    
    #Create first line of sensibility matrix
    BTTB = bttb(x,y,z,h)
    
    #Calculates the eigenvalues of BCCB matrix
    cev = bccb(shape,N,BTTB)

    #Fast Equivalent layer loop
    m_new = fast_loop_bccb(cev,shape,N,data,diagonal_A,rho0,itmax)

    #Final predicted data
    gzp = fast_forward_bccb(shape,N,m_new,cev)

    return m_new, gzp

def diagonal(x,y,shape):
    '''
    Calculates a NxN diagonal matrix given by
    the area of x and y data spacing.

    input
    x, y: numpy array - the x, y coordinates of
    the grid and equivalent layer points.
    shape: tuple - grid size.

    output
    N: scalar - number of observation points.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.
    '''
    N = shape[0]*shape[1]
    assert N == x.size, 'N and x must have the same number\
    of elements'
    delta_s = ((np.max(x)-np.min(x))/(shape[0]-1))*\
    ((np.max(y)-np.min(y))/(shape[1]-1))
    diagonal_A = G*SI2MGAL*2.*np.pi/(delta_s)
    return N,diagonal_A

def diagonal_gzz(x,y,shape):
    '''
    Calculates a NxN diagonal matrix given by
    the area of x and y data spacing.

    input
    x, y: numpy array - the x, y coordinates of
    the grid and equivalent layer points.
    shape: tuple - grid size.

    output
    N: scalar - number of observation points.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.
    '''
    N = shape[0]*shape[1]
    assert N == x.size, 'N and x must have the same number\
    of elements'
    delta_s = ((np.max(x)-np.min(x))/(shape[0]-1))*\
    ((np.max(y)-np.min(y))/(shape[1]-1))
    diagonal_A = G*10**(6)*2*np.pi/(delta_s)
    return N,diagonal_A

def i_rho(data,diagonal_A):
    '''
    Calculates the initial equivalent layer
    property estimative.

    input
    data: numpy array - the gravity disturbance.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.

    output
    rho0: numpy array - initial equivalent
    layer property estimative.
    '''
    rho0 = data/diagonal_A
    return rho0

def bttb(x,y,z,h):
    '''
    Calculates the first line of sensbility matrix.

    input
    x, y: numpy array - the x, y coordinates
    of the grid and equivalent layer points.
    z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.

    output
    W_bt: numpy array - first line os sensibility matrix.
    '''
    a = (x-x[0])
    b = (y-y[0])
    c = (h-z[0])
    W_bt = (G*SI2MGAL*c)/((a*a+b*b+c*c)**(1.5))
    return W_bt

def bccb(shape,N,BTTB):
    '''
    Calculates the eigenvalues of the BCCB matrix.

    input
    shape: tuple - grid size.
    N: scalar - number of observation points.
    BTTB: numpy array - first line os sensibility matrix.

    output
    cev: numpy array - eigenvalues of the BCCB matrix.
    '''
    cev = np.zeros(4*N, dtype='complex128')
    k = 2*shape[0]-1
    for i in xrange (shape[0]):
        block = BTTB[shape[1]*(i):shape[1]*(i+1)]
        rev = block[::-1]
        cev[shape[1]*(2*i):shape[1]*(2*i+2)] = np.concatenate((block,0,rev[:-1]), axis=None)
        if i > 0:
            cev[shape[1]*(2*k):shape[1]*(2*k+2)] = cev[shape[1]*(2*i):shape[1]*(2*i+2)]
            k -= 1
    cev = cev.reshape(2*shape[0],2*shape[1]).T
    return np.fft.fft2(cev)

def bttb_gzz(x,y,z,h):
    '''
    Calculates the first line of sensbility matrix.

    input
    x, y: numpy array - the x, y coordinates
    of the grid and equivalent layer points.
    z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.

    output
    W_bt: numpy array - first line os sensibility matrix.
    '''
    a = (x-x[0])
    b = (y-y[0])
    c = (h-z[0])
    W_bt = -1/((a*a+b*b+c*c)**(1.5)) + (3*c*c)/((a*a+b*b+c*c)**(2.5))
    W_bt = G*10**(9)*W_bt
    return W_bt

def sensibility_matrix(x,y,z,h,N):
    '''
    Calculates a full NxN matrix given by
    the first derivative of the function
    1/r.

    input
    x, y: numpy array - the x, y coordinates of
    the grid and equivalent layer points.
	z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.
    N: scalar - number of observation points.

    output
    A: matrix - full NxN matrix given by
    the first derivative of the function
    1/r.
    '''
    A = np.zeros((N, N), dtype=np.float)

    for i in xrange (N):
        a = (x-x[i])
        b = (y-y[i])
        c = (h-z[i])
        A[i] = c/((a*a+b*b+c*c)**(1.5))
    A = A*G*SI2MGAL
    return A

def sensibility_matrix_gzz(x,y,z,h,N):
    '''
    Calculates a full NxN matrix given by
    the first derivative of the function
    1/r.

    input
    x, y: numpy array - the x, y coordinates of
    the grid and equivalent layer points.
	z: numpy array - the height of observation points.
    h: numpy array - the depth of the equivalent layer.
    N: scalar - number of observation points.

    output
    A: matrix - full NxN matrix given by
    the first derivative of the function
    1/r.
    '''
    A = np.zeros((N, N), dtype=np.float)

    for i in xrange (N):
        a = (x-x[i])
        b = (y-y[i])
        c = (h-z[i])
        A[i] =  -1/((a*a+b*b+c*c)**(1.5)) + 3*c*c/((a*a+b*b+c*c)**(2.5))
    A = A*G*SI2MGAL
    return A

def fast_loop(data,A,diagonal_A,rho0,itmax):
    '''
    Solves the linear inversion through a iterative method.

    input
    data: numpy array - the gravity disturbance.
	A: matrix - full NxN matrix given by
    the first derivative of the function
    1/r.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.
	rho0: numpy array - initial equivalent
    layer property estimative.
	itmax: scalar - number of iterations

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    '''
    m_new = np.copy(rho0)
    for i in xrange (itmax):
        res = (data - A.dot(m_new))
        delta_m = res/diagonal_A
        m_new += delta_m
    return m_new

def fast_loop_toep(BTTB,shape,N,data,diagonal_A,rho0,itmax):
    '''
    Solves the linear inversion through a iterative method.

    input
	BTTB: numpy array - first line os sensibility matrix.
	shape: tuple - grid size.
	N: scalar - number of observation points.
    data: numpy array - the gravity disturbance.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.
	rho0: numpy array - initial equivalent
    layer property estimative.
	itmax: scalar - number of iterations

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    '''
    m_new = np.copy(rho0)
    for i in xrange (itmax):
        gzp = fast_forward_toep(shape,N,m_new,BTTB)
        res = (data - gzp)
        delta_m = res/diagonal_A
        m_new += delta_m
    return m_new
    
def fast_loop_bccb(cev,shape,N,data,diagonal_A,rho0,itmax):
    '''
    Solves the linear inversion through a iterative method.

    input
	BTTB: numpy array - first line os sensibility matrix.
	shape: tuple - grid size.
	N: scalar - number of observation points.
    data: numpy array - the gravity disturbance.
    diagonal_A: escalar - NxN diagonal matrix
    given by the area of x nd y data spacing.

    output
    m_new: numpy array - final equivalent
    layer property estimative.
    '''
    m_new = np.copy(rho0)
    for i in xrange (itmax):
        gzp = fast_forward_bccb(shape,N,m_new,cev)
        res = (data - gzp)
        delta_m = res/diagonal_A
        m_new += delta_m
    return m_new

def fast_forward_toep(shape,N,p,BTTB):
    '''
    Calculate the forward problem at each iteration
    taking advantage of the BTTB (Block-Toeplitz
    Toreplitz-Block) structures.

    input
	shape: tuple - grid size.
    N: scalar - number of observation points.
    p: numpy array - equivalent layer property
	estimative.
	BTTB: numpy array - first line os sensibility matrix.

    output
    gzp: numpy array - the predicted data.
    '''
    gzp = np.zeros(N).reshape(N,1)
    p = p.reshape(N,1)

    l = shape[0]

    for i in xrange (shape[0]):
        k = 0 #Index of the parameter's segment
        j = 0+i #Index of the data's segment using the lower matrix
        u = 0 #Index of the data's segment using the upper matrix

        #Create each Toeplitz block by the segment of the first row
        block = toeplitz(BTTB[shape[1]*(i):shape[1]*(i+1)])

        while k < l: # lower matrix
            block_p = block.dot(p[shape[1]*(k):shape[1]*(k+1)])
            gzp[shape[1]*(j):shape[1]*(j+1)] += block_p
            if k >= i and i > 0: # upper matrix 1
                gzp[shape[1]*(u):shape[1]*(u+1)] += block_p
                u += 1
            k += 1
            j += 1
        while i > 0 and k >= l and k < shape[0]: # upper matrix 2
            if k >= i:
                block_p = block.dot(p[shape[1]*(k):shape[1]*(k+1)])
                gzp[shape[1]*(u):shape[1]*(u+1)] += block_p
                u += 1
            k += 1
        l -= 1

    return np.ravel(gzp)
    
def fast_forward_bccb(shape,N,p,cev):
    '''
	Calculate the forward problem at each iteration
    taking advantage of the BTTB (Block-Toeplitz
    Toreplitz-Block) structures.

    input
    shape: tuple - grid size.
    N: scalar - number of observation points.
    p: numpy array - equivalent layer property
	estimative.
	BTTB: numpy array - first line os sensibility matrix.

    output
    gzp: numpy array - the predicted data.
    '''
    v = np.zeros(4*N, dtype='complex128')
    for i in xrange (shape[0]):
        v[shape[1]*(2*i):shape[1]*(2*i+2)] = np.concatenate((p[shape[1]*(i):shape[1]*(i+1)], np.zeros(shape[1])), axis=None)
    
    v = v.reshape(2*shape[0],2*shape[1]).T
    gzp = np.fft.ifft2(np.fft.fft2(v)*cev)
    gzp = np.ravel(np.real(gzp[:shape[1],:shape[0]]).T)
    return gzp

def classic_mag(x,y,z,zj,F,h,N,data):
    '''
    '''
    A = np.empty((N, N), dtype=np.float)
    for i in xrange (N):
        a = (x-x[i])
        b = (y-y[i])
        c = (zj-z[i])
        r = (a*a+b*b+c*c)
        r3 = r**(-1.5)
        r5 = r**(2.5)
        Hxx = -r3+3*(a*a)/r5
        Hxy = 3*(a*b)/r5
        Hxz = 3*(a*c)/r5
        Hyy = -r3+3*(b*b)/r5
        Hyz = 3*(b*c)/r5
        Hzz = -r3+3*(c*c)/r5
        A[i] = 100*((F[0]*Hxx+F[1]*Hxy+F[2]*Hxz)*h[0] + (F[0]*Hxy+F[1]*Hyy+F[2]*Hyz)*h[1] + (F[0]*Hxz+F[1]*Hyz+F[2]*Hzz)*h[2])
    I = np.identity(N)
    ATA = A.T.dot(A)
    mu = (np.trace(ATA)/N)*10**(-2)
    AI = inv(ATA+mu*I)
    p = (AI.dot(A.T)).dot(data)
    tf = A.dot(p)
    return p, tf 
    
def bttb_mag(x,y,z,zj,F,h,shape):
    '''
    '''
    # Fisrt line of BTTB first row
    a = (x-x[0])
    b = (y-y[0])
    c = (zj-z[0])
    r = (a*a+b*b+c*c)
    r3 = r**(-1.5)
    r5 = r**(2.5)
    Hxx = -r3+3*(a*a)/r5
    Hxy = 3*(a*b)/r5
    Hxz = 3*(a*c)/r5
    Hyy = -r3+3*(b*b)/r5
    Hyz = 3*(b*c)/r5
    Hzz = -r3+3*(c*c)/r5
    bttb_0 = 100*((F[0]*Hxx+F[1]*Hxy+F[2]*Hxz)*h[0] + (F[0]*Hxy+F[1]*Hyy+F[2]*Hyz)*h[1] + (F[0]*Hxz+F[1]*Hyz+F[2]*Hzz)*h[2])
    
    # Last line of BTTB first row
    a = (x-x[shape[1]-1])
    b = (y-y[shape[1]-1])
    c = (zj-z[shape[1]-1])
    r = (a*a+b*b+c*c)
    r3 = r**(-1.5)
    r5 = r**(2.5)
    Hxx = -r3+3*(a*a)/r5
    Hxy = 3*(a*b)/r5
    Hxz = 3*(a*c)/r5
    Hyy = -r3+3*(b*b)/r5
    Hyz = 3*(b*c)/r5
    Hzz = -r3+3*(c*c)/r5
    bttb_1 = 100*((F[0]*Hxx+F[1]*Hxy+F[2]*Hxz)*h[0] + (F[0]*Hxy+F[1]*Hyy+F[2]*Hyz)*h[1] + (F[0]*Hxz+F[1]*Hyz+F[2]*Hzz)*h[2])
    
    # Fisrt line of BTTB last row
    a = (x-x[-shape[1]])
    b = (y-y[-shape[1]])
    c = (zj-z[-shape[1]])
    r = (a*a+b*b+c*c)
    r3 = r**(-1.5)
    r5 = r**(2.5)
    Hxx = -r3+3*(a*a)/r5
    Hxy = 3*(a*b)/r5
    Hxz = 3*(a*c)/r5
    Hyy = -r3+3*(b*b)/r5
    Hyz = 3*(b*c)/r5
    Hzz = -r3+3*(c*c)/r5
    bttb_2 = 100*((F[0]*Hxx+F[1]*Hxy+F[2]*Hxz)*h[0] + (F[0]*Hxy+F[1]*Hyy+F[2]*Hyz)*h[1] + (F[0]*Hxz+F[1]*Hyz+F[2]*Hzz)*h[2])
    
    # Last line of BTTB last row
    a = (x-x[-1])
    b = (y-y[-1])
    c = (zj-z[-1])
    r = (a*a+b*b+c*c)
    r3 = r**(-1.5)
    r5 = r**(2.5)
    Hxx = -r3+3*(a*a)/r5
    Hxy = 3*(a*b)/r5
    Hxz = 3*(a*c)/r5
    Hyy = -r3+3*(b*b)/r5
    Hyz = 3*(b*c)/r5
    Hzz = -r3+3*(c*c)/r5
    bttb_3 = 100*((F[0]*Hxx+F[1]*Hxy+F[2]*Hxz)*h[0] + (F[0]*Hxy+F[1]*Hyy+F[2]*Hyz)*h[1] + (F[0]*Hxz+F[1]*Hyz+F[2]*Hzz)*h[2])
    
    return bttb_0, bttb_1, bttb_2, bttb_3

def fast_forward_bccb_mag(bttb_0,bttb_1,bttb_2,bttb_3,p,shape,N):
    '''
    '''
    # First column of BCCB
    BCCB = np.zeros(4*N, dtype='complex128')
    v = np.zeros(4*N, dtype='complex128')
    q = shape[0]-1
    k = 2*shape[0]-1
    for i in xrange (shape[0]):
        v[shape[1]*(2*i):shape[1]*(2*i+2)] = np.concatenate((p[shape[1]*(i):shape[1]*(i+1)], np.zeros(shape[1])), axis=None)
        
        block_2 = bttb_2[shape[1]*(q):shape[1]*(q+1)]
        block_3 = bttb_3[shape[1]*(q):shape[1]*(q+1)]
        c_2 = block_2[::-1]
        BCCB[shape[1]*(2*i):shape[1]*(2*i+2)] = np.concatenate((block_3[::-1],0,c_2[:-1]), axis=None)
        q -= 1
    
        if i > 0:
            block_0 = bttb_0[shape[1]*(i):shape[1]*(i+1)]
            block_1 = bttb_1[shape[1]*(i):shape[1]*(i+1)]
            c_0 = block_0[::-1]
            BCCB[shape[1]*(2*k):shape[1]*(2*k+2)] = np.concatenate((block_1[::-1],0,c_0[:-1]), axis=None)
            k -= 1
    
    BCCB = BCCB.reshape(2*shape[0],2*shape[1]).T
    v = v.reshape(2*shape[0],2*shape[1]).T
    #Matrx-vector product
    dobs_bccb = np.fft.ifft2(np.fft.fft2(v)*np.fft.fft2(BCCB))
    dobs_bccb = np.ravel(np.real(dobs_bccb[:shape[1],:shape[0]]).T)
    return dobs_bccb