%%
N = 4;

% Forcing to be odd
if ( mod(N+1, 2) )  
    N = N + 1;
end

%% Points %
t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );

[rr, zz, tt] = ndgrid(r, z, t);

disp("printing rr");
disp(rr);
disp("printing zz");
disp(zz);
disp("printing tt");
disp(tt);

%% Matrices

Ops = Operators(N,'C02, D2, A, coF');
C02 = Ops.C02;
D2 = Ops.D2;
A = Ops.A;
coF = Ops.coF;

%% Initial Condition
CFS = zeros(size(zz));
f = @(r,z,th,t) r .* sin(th) + z .* cos(th);
% finit = @(r,z,th,t) f(r,z,th,0);


for j=1:N
    CFS(:,j,:) = f(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

disp("printing CFS:");
disp(CFS);