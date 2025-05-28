function X = func2grid(f,N,ord,coeff)
% This code takes a function of K variables as its argument, and samples it
% at a grid of size NxNx...xN. The shape of the grid is defined by the
% string 'ord': the character 't' denotes a trigonometric dimension, the 
% character 'e' denotes an equispaced dimension (from -1 to 1), and any 
% other letter denotes a Chebyshev dimension. For instance, if ord = 'rzt', 
% then the grid is Chebyshev by Chebyshev by trigonometric. If coeff is
% given (and not 0), the output grid will be transformed to trigonometric
% and Chebyshev coefficient space. This will fail for equispaced grids.

% Dave Darrow. February 3, 2018.

t = pi*trigpts( N );
r = chebpts( N );
eq = linspace(-1,1,N).';

grid = cell(length(ord),1);
for i = 1:length(ord)
    if ( ord(i) == 't' || ord(i) == 'T' )
        grid{i} = t;
    elseif ( ord(i) == 'e' || ord(i) == 'E' )
        grid{i} = eq;
    else
        grid{i} = r;
    end
end
GRID = cell(length(ord),1);
[GRID{:}] = ndgrid(grid{:});
X = f(GRID{:});

if ( nargin > 3 && coeff )
    X = V2C_cyl(X,ord);
end
end