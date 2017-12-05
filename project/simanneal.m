function [x_min, fx_min] = ...
    simanneal( f, a, b, cooling_factor, nrand, maxiter, min_gamma, ...
               max_non_improving_steps, alpha, verbose )
  %--------------------------------------------------------------------------
  % Jonathan R. Senning <jonathan.senning@gordon.edu>
  % Gordon College
  % April 20, 2005
  % Revised November 26, 2008 Added cf (cooling factor) parameter
  % Revised December  9, 2008 Added nr (number of random trials) parameter
  % Revised December 11, 2008 Adapted to work with both Octave and Matlab
  %
  % Usage: [x_min, fx_min] = simanneal(@f, a, b, ...)
  %
  % Returns:
  %       x_min:     Location of minimum.
  %
  %       fx_min:    Minimum value of function.
  %
  % Parameters:
  %       f:         Handle of the objective function to minimize.
  %
  %       a:         The left endpoint of the interval.
  %
  %       b:         Right endpoint of interval.
  %
  %       cf:        Cooling factor
  %
  %       nrand:     Number of random values per trial
  %
  %       maxiter:   Maximum number of iterations to carry out
  %
  %       min_gamma: Minimum width of random distribution
  %
  %       max_non_improving_steps: Maximum number of non-improving steps
  %                  to carry out before gamma is reduced
  %
  %       alpha:     Likelihood that a higher uphill move is made in
  %                  preference to a lower uphill move.
  %
  %       verbose:   If nonzero then information about each step is displayed.
  %
  % Uses simulated annealing to find the global minimum of a function of one
  % variable.  This was written sometime in the past; unfortunately I don't
  % remember where I got the algorithm from :(.
  %--------------------------------------------------------------------------

  % Define default values for optional arguments.  The order here must
  % reflect the order they appear in the parameter list.

  if ( nargin < 10 )
    % Nonzero (true) value indicates the user wants to see output from
    % each iteration.
    verbose = 0;
  end
  if ( nargin < 9 )
    % This positive parameter determines the likelihood that a higher
    % uphill move is made over a lower uphill move.
    alpha = 1.0;
  end
  if ( nargin < 8 )
    % The number of iterations that fail to produce an improvement before
    % the system will be cooled (i.e. the value of gamma is reduced by a
    % factor of cooling_factor).
    max_non_improving_steps = 5;
  end
  if ( nargin < 7 )
    % Smallest value of gamma before iteration terminates.
    min_gamma = 1e-8;
  end
  if ( nargin < 6 )
    % The main iteration terminates when gamma becomes smaller than
    % min_gamma or when the number of iterations exceeds maxiter.
    maxiter = 10000;
  end
  if ( nargin < 5 )
    % This is the maximum number of random trials performed during each
    % iteration; the rejection of values outside the interval may
    % eliminate some.  This means the objective function is evaluated at
    % most this many times during each iteration.
    nrand = 10;
  end
  if ( nargin < 4 )
    % This should be in the range (0,1).  Values close to 1 result in slow
    % cooling and should allow the best chance of finding the true minimum
    % (analogous to the minimum energy state) but will result in slow
    % convergence.  Values close to 0 will result in fast convergence to a
    % local minimum, but likely not the absolute minimum.
    cooling_factor = 0.8;
  end

  % Initialize for the iterations.  Gamma controls the width of the normally
  % distributed random numbers used in each iteration.  By initializing to
  % (b-a)/2 we start with random numbers distributed over the entire interval.

  gamma = 0.5 * ( b - a );
  x  = 0.5 * ( a + b );
  fx = f( x );
  x_min  = x;
  fx_min = fx;
  non_improving_steps = 0;

  % Main loop

  for k = 1 : maxiter

    % Generate normally distributed random number centered around x and then
    % reject any that lie outside the interval [a,b].

    u = sort( gamma * randn( nrand, 1 ) + x );
    lu = length( u );
    m = 1;
    while ( m <= lu && u(m) < a )
      m = m + 1;
    end
    n = nrand;
    while ( n >= 1 && u(n) > b )
      n = n - 1;
    end
    nr = n - m + 1;
    u = u(m:n);
    if ( length( u ) == 0 )
      % no random numbers remain after trimming, skip this iteration
      continue
    end

    % Determine the function value at each of the random points, find the
    % minimum value and also the random value that gave it.  The second
    % parameter returned by min() is the index of the minimum value in the
    % fu array.

    fu = f( u );
    [fu_min, i] = min( fu );
    u_min = u(i);

    % Have we found an improvement over our current center point?  If so, use
    % it; if not, select a new center point chosen from the current random
    % distribution.

    if ( fu_min < fx )
      x  = u_min;
      fx = fu_min;
    else
      p = exp( alpha * ( fx - fu ) );
      p = p / sum( p );
      t = rand;
      i = 1;
      s = p(i);
      while ( s < t && i < nr )
        i = i + 1;
        s = s + p(i);
      end
      x  = u(i);
      fx = fu(i);
    end

    % See if we've found an improved minimum.  If we've not found one in a
    % while, reduce gamma (i.e. cool the system) and try again.

    if ( fx < fx_min )
      x_min  = x;
      fx_min = fx;
      non_improving_steps = 0;
    else
      non_improving_steps = non_improving_steps + 1;
      if ( non_improving_steps > max_non_improving_steps )
        non_improving_steps = 0;
        gamma = cooling_factor * gamma;
        x  = x_min;
        fx = fx_min;
      end
    end

    % Done with iteration.  Display status (if requested) and check
    % termination criteria.

    if ( verbose ~= 0 )
      fprintf( '%4d: %20.16f %20.16f %8.5f %8.5f  %2d %4.2e\n', ...
               k, x_min, fx_min, x, fx, non_improving_steps, gamma );
      %fflush( 1 );
    end

    if ( gamma < min_gamma )
      break;
    end

  end

% End of file
