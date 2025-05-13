%%
N = 4;
T = 1;

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
finit = @(r,z,th) f(r,z,th,0);


for j=1:N
    CFS(:,j,:) = finit(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

disp("printing CFS:");
disp(CFS);

%% V2C

% from trigtech and chebtech2
CFS_func = V2C_cyl(CFS,'rzt');
CFS(:,:,:,2:T)=zeros; % sets all T > 1 to be 0, when T=1 it does nothing
disp('CFS_func: ');
disp(CFS_func); % prints N many columns per every T

% written explicitly
CFS_expl = V2C_cyl_manual(CFS,'rzt');
CFS(:,:,:,2:T)=zeros;
disp('CFS_expl: ');
disp(CFS_expl .* -1);

%% Comparing
vals_cheb = [-1; -0.5; 0; 0.5; 1];
display(vals_cheb)
vals_cheb_updated = chebv2c(vals_cheb);
disp("vals updated")
display(vals_cheb_updated)

vals_trig = [-1, -0.5i, 0, 0.5i, 1; -1i, -0.5, 0, 0.5, 1i; -1, -0.5, 0i, 0.5, 1];
display(vals_trig);
vals_updated = trigv2c(vals_trig);
disp("vals updated")
display(vals_updated);

%% Testing Trigtech Function

vals = pi* linspace(-1, 1, N);
disp('vals:');
disp(vals);

% this is the chebfun function!
ttfunc = trigtech.vals2coeffs(vals);

% this is the explicit equation!
function c = vals2coeffs(vals)
    n = length(vals);               % Number of sample points
    c = fftshift(fft(vals))/n;       % Compute FFT and normalize
end


ttexpl = vals2coeffs(vals);

disp('ttfunc:');
disp(ttfunc);

disp('ttexpl:');
disp(ttexpl);


