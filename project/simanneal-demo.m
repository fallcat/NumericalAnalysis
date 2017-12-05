#!/usr/bin/octave -q

%-----------------------------------------------------------------------------
% Jonathan R. Senning <jonathan.senning@gordon.edu>
% Gordon College
% April 20, 2005
% Revised April 27, 2006 to update graphing commands
% Revised November 26, 2008 to update commands to save PostScript file
% Revised December 11, 2008 to work with both Octave and MatLab
%
% $Id: simanneal-demo,v 1.5 2007/03/26 20:02:12 senning Exp senning $
%
% Demonstrate how simulated annealing is used to find the global minimum of a
% function with multiple local minima.  The function simanneal() was written
% sometime in the past; unfortunately I don't remember where I got the
% algorithm from :(.
%-----------------------------------------------------------------------------

% Define interval and function

a = 0;
b = 5;
cooling_factor = 0.8;
f = @(x) x .* exp(-0.5*x) .* sin(55*x) .* cos(12*x);

% Perform simulated annealing.  The last parameter may be changed to a
% nonzero value to print out intermediate results while searching for the
% minimum.

[x_min, fx_min] = simanneal( f, a, b, cooling_factor );
fprintf( 'f(%f) = %f\n', x_min, fx_min );

% Display graph of results.  The graph of the function f(x) is shown on the
% interval [a,b] and the minimum point is plotted as a single point.

x = linspace( a, b, 1001 );
y = f( x );

plot( x, y, 'b-', x_min, fx_min, 'ro' );
xlabel( 'x' );
ylabel( 'y' );
title( 'Minimization of Continuous Function by Simulated Annealing' );
legend( 'Function', 'Minimum', 'location', 'NorthEast' );
drawnow();

yesno='';
yesno = input( 'Create PNG file? [y/N] ', 's' );

if ( length( yesno ) > 0 && ( yesno(1) == 'Y' || yesno(1) == 'y' ) )
  pngfile = 'simanneal.png';
  print( pngfile, '-dpng' );
  fprintf( 'PNG file %s was created\n', pngfile );
end

% End of File
