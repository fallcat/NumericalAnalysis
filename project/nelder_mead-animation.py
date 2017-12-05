#!/usr/bin/env python

import pylab
import numpy as np
import time

def draw_domain( f, x0, y0, x1, y1 ):
    number_of_levels = 41
    n = 101
    x = pylab.linspace( x0, x1, n )
    y = pylab.linspace( y0, y1, n )
    xx = ( pylab.reshape(x,(n,1)) * pylab.ones(n) ).transpose()
    yy = pylab.reshape(y,(n,1)) * pylab.ones(n)
    z = f( [xx, yy] )
    pylab.ion()
    pylab.contour( x, y, z, number_of_levels )
    pylab.draw()

def nelder_mead( f, x, dx, maxerr, maxiter, verbose = False ):
    """Finds minimum of multivarite function f(x)

    USAGE:
        y, err, iter = nelder_mead( f, x, dx, maxerr, maxiter, verbose = False )

    INPUT:
        f       - The objective function to minimize; must accept a single
                  argument which is an n-element array.
        x       - Initial estimate of the minimum point; should be supplied
                  as an n-element list or array.
        dx      - Displacement vector for the vertices in the initial
                  simplex.  All but one vertex of the simplex are displaced
                  in one dimension from the initial point using the values
                  in this array.  May either be a scalar value, in which case
                  it applies to each dimension of the simplex, or an
                  n-element list or array.  A reasonable choice is 0.1 * x.
        maxerr  - Maximum allowable relative difference between the best and
                  worst vertices in the simplex.  This provides the stopping
                  criteria and may either be a scalar value, in which case it
                  applies to each vertex of the simplex, or an n-element
                  vector.
        maxiter - Maximum number of iterations to perform.  Procedure will
                  terminate after this many iterations are done regardless
                  of the error.
        verbose - Optional Parameter -- if True then information about each
                  step is displayed.

    OUTPUT:
        y       - Array containing coordinates of point which minimizes the
                  objective function.
        err     - Array containing the relative difference between the the
                  best estimate and the worst estimate of the minimum point
                  in the final simplex.
        iter    - Number of iterations performed.

    NOTES:
        This function uses the Simplex method of Nelder and Mead to search
        for the minimum of a function f(x) where x is an n-element row
        vector that contains the values of the n independent variables used
        by f.

        This code is based, through several derivations, on code presented in
        "Fitting Curves to Data", by Marco Caceci and William Cacheris, Byte
        Magazine, May 1984.

    AUTHOR:
        Jonathan R. Senning <jonathan.senning@gordon.edu>
        Gordon College
        Original Pascal Version: 1985
        C Version: 1999
        Octave Version: May 7, 2003
        Python Version: Mar 13, 2008
    """

    # define algorithm parameters

    alpha = 2.0  # reflection factor
    beta =  0.5  # contraction factor
    gamma = 3.0  # expansion factor
    delta = 0.5  # shrinkage factor
    epsilon = np.finfo( np.float ).eps  # machine epsilon

    # find number of values we are searching for make sure they are in
    # a NumPy array

    n = len( x )
    x = np.array( x, float )

    # convert dx and maxerr to NumPy arrays if they are passed as scalars

    try:
        size = len( dx )
    except TypeError:
        dx = [dx] * n
    dx = np.array( dx, float )

    try:
        size = len( maxerr )
    except TypeError:
        maxerr = [maxerr] * n
    maxerr = np.array( maxerr, float )
    err    = 2.0 * maxerr

    # construct initial simplex.  each row represents a vertex in
    # (n+1)-dimensional space and we need (n+1) vertices.

    m = n + 1
    simplex = np.array( [x] * m, float )
    value = np.zeros( m )

    # adjust the first n vertices using the dx offsets and compute the
    # value of the objective function at each simplex vertex

    for i in xrange( n ):
        simplex[i,i] += dx[i]

    for i in xrange( m ):
        value[i] = f( simplex[i] )

    # identify best (smallest objective function value) and worst (largest
    # objective function value) vertices

    best, worst = ( value.argmin(), value.argmax() )

    if verbose:
        print "Initial Simplex"
        print np.append( simplex, value.reshape( m, 1 ), 1 )
        ux = [0, 0];
        vx = [0, 0];

    # start main loop.  iterations will continue as long as the relative
    # difference between best and worst objective values in the simplex
    # exceeds the maximum allowable error for each vertex or until the
    # maximum number of iterations is reached.

    iter = 0
    while ( err > maxerr ).any() and iter < maxiter:

        iter = iter + 1

        if verbose: print "Iteration %5d: " % iter,

        # compute centroid of all vertices in the simplex but the worst
        center = sum( simplex[[i for i in range( m ) if i != worst]] ) / n

        # try reflection
        u = alpha * center + ( 1 - alpha ) * simplex[worst]
        reflection = f( u )

        if reflection <= value[best]:

            # reflection was good - perhaps we can do better: try expansion
            w = gamma * center + ( 1 - gamma ) * simplex[worst]
            expansion = f( w )

            if expansion <= value[best]:

                # expansion is good, keep it
                simplex[worst] = w
                value[worst] = expansion
                if verbose: print "Expansion"

            else:

                # expansion not good enough, keep reflection
                simplex[worst] = u
                value[worst] = reflection
                if verbose: print "Reflection"

        else:

            # compare reflection to worst candidate
            if reflection <= value[worst]:

                # reflection is an improvement, keep it
                simplex[worst] = u
                value[worst] = reflection
                if verbose: print "Reflection"

            else:

                # reflection is worse, try contraction
                w = beta * center + ( 1 - beta ) * simplex[worst]
                contraction = f( w )

                if contraction <= value[worst]:

                    # contraction is good enough to keep
                    simplex[worst] = w
                    value[worst] = contraction
                    if verbose: print "Contraction"

                else:

                    # contraction not good enough
                    # shrink all non-best vertices toward best vertex
                    if verbose: print "Shrinking"
                    for i in [j for j in xrange( m ) if j != best]:
                        simplex[i] = ( 1 - delta ) * simplex[i] \
                            + delta * simplex[best]
                        value[i] = f( simplex[i] )

        # done with adjustments, update best and worst vertices

        best, worst = ( value.argmin(), value.argmax() )

        if verbose:
            pylab.plot( ux, vx, 'k-' )
            print np.append( simplex, value.reshape( m, 1 ), 1 )
            ux = np.array( simplex )[:,0]
            ux = np.append( ux, ux[0] )
            vx = np.array( simplex )[:,1]
            vx = np.append( vx, vx[0] )
            print ux, vx
            pylab.plot( ux, vx, 'r-' )
            pylab.draw()
            time.sleep( 0.5 )

        # compute relative error at each vertex; make sure we don't divide
        # by zero

        denom = simplex[best]
        abs_denom = abs( denom )
        while abs_denom.min() == 0.0:
            denom[abs_denom.argmin()] = epsilon
        err = abs( ( simplex[worst] - simplex[best] ) / denom )

        # end of main while loop

    # return the best vertex as our result

    return ( simplex[best], err, iter )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Example of how the nelder_mead() function can be used.  In this case it
# is used to perform nonlinear least-squares curvefitting.

if __name__ == "__main__":

    def f( x ):
	return (x[0]-1.0)**2 + (x[1]-2.0)**2 + 0.8 * np.cos( x[0] * x[1] )


    n = 101
    x = pylab.linspace( -2.0, 4.0, n )
    y = pylab.linspace( -1.0, 5.0, n )
    xx = ( pylab.reshape(x,(n,1)) * pylab.ones(n) ).transpose()
    yy = pylab.reshape(y,(n,1)) * pylab.ones(n)
    z = f( [xx, yy] )

#    x0 = [3.8, 3.5]
#    dx = 0.25
#    x0 = [3.8, 0.0]
#    dx = 0.75
    x0 = [-1.5, 3.5]
    dx = 0.25

    pylab.ion()
    pylab.contour( x, y, z, 41 )
    a, err, iter = nelder_mead( f, x0, dx, 1e-3, 1, True )
    pylab.draw()

    raw_input( 'Press Enter to start...' )

    # search for minimum of ssqr() function

    a, err, iter = nelder_mead( f, x0, dx, 1e-2, 100, True )

    # report results

    print "a    = ", a
    print "err  = ", err
    print "iter = ", iter

    raw_input( 'Press Enter to quit...' )
