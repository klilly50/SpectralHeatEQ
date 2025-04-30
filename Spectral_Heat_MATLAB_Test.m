
%% 

N = 10;
% Forcing to be even

if ( mod(N+1, 2) )  
    N = N + 1;
end

%% Points %

N=4;

t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );

[rr, zz, tt] = ndgrid(r, z, t);
