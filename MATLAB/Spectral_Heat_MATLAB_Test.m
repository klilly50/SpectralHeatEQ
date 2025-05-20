%%
N = 4;
T = 10;

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

%% V2C ORIGINAL
disp('FUNCTION')
CFS = V2C_cyl(CFS,'rzt');
CFS(:,:,:,2:T)=zeros; % sets all T > 1 to be 0, when T=1 it does nothing

disp('CFS: ');
disp(CFS);
disp('CFS size: ');
disp(size(CFS));

% %% V2C EXPLICILY WRITTEN
% disp('EXPLICIT')
% CFS_expl = V2C_cyl_manual(CFS,'rzt');
% CFS(:,:,:,2:T)=zeros;
% 
% disp('CFS_expl: ');
% disp(CFS_expl);
% % disp(CFS_expl .* -1);
% 
% %% Comparing
% vals_cheb = [-1; -0.5; 0; 0.5; 1];
% display(vals_cheb)
% vals_cheb_updated = chebv2c(vals_cheb);
% disp("vals updated")
% display(vals_cheb_updated)
% 
% vals_trig = [-1, -0.5i, 0, 0.5i, 1; -1i, -0.5, 0, 0.5, 1i; -1, -0.5, 0i, 0.5, 1];
% display(vals_trig);
% vals_updated = trigv2c(vals_trig);
% disp("vals updated")
% display(vals_updated);
% 
% %% Testing Trigtech Function
% 
% vals = pi* linspace(-1, 1, N);
% disp('vals:');
% disp(vals);
% 
% % this is the chebfun function!
% ttfunc = trigtech.vals2coeffs(vals);
% 
% % this is the explicit equation!
% function c = vals2coeffs(vals)
%     n = length(vals);               % Number of sample points
%     c = fftshift(fft(vals))/n;       % Compute FFT and normalize
% end
% 
% 
% ttexpl = vals2coeffs(vals);
% 
% disp('ttfunc:');
% disp(ttfunc);
% 
% disp('ttexpl:');
% disp(ttexpl);

%% Plug in first four time steps
alph = 1e-9;
Hist = cell(3,1);

for j = 1:3
    for k=1:N
        Hist{j}(:,k,:) = f(rr(:,k,:),zz(:,k,:), tt(:,k,:),j*alph);
    end
    CFS(:,:,:,j+1)=V2C_cyl(Hist{j},'rzt');
end

disp(CFS)

%% Full Time Step
q = (12/25)*alph;
for d = 5:T
        RHS = (48/25)*CFS(:,:,:,d-1)-(36/25)*CFS(:,:,:,d-2)+(16/25)*CFS(:,:,:,d-3)-(3/25)*CFS(:,:,:,d-4);
        CFS(:,:,:,d) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

% Convert back to value space
% CFS(:,:,:,T) = C2V_cyl(CFS(:,:,:,T),'rzt');

%% Screened Poisson Solver

function Fapprox = Heat3DSolver(cfs,Q,N,C02,A,D2,coF)

    % Solves the screened Poisson equation, [I-Q*Del^2]u=CFS in coeff. space
    % Homogeneous Dirichlet boundary conditions.
    
    % Discretization size: N
    
    Fapprox = zeros(N,N,N);
    
    bc = [ (-1).^(0:N-1) > 0 ; 
           (-1).^(0:N-1) < 0 ]; 
    
    for k = -(N-1)/2 : (N-1)/2
        
        AA = coF - Q*(A - k^2*C02);
        BB = C02; 
        CC = -Q*coF; 
        DD = D2; 
        RHS = coF*cfs(:,:,k+(N+1)/2)*C02.';
        %RHS = RHS(1:N-2, 1:N-2);
        
        [AA, RHS] = zeroDOF(AA, BB, RHS, bc, zeros(2,N));
        [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, zeros(2,N));
        RHS = RHS.';
        [CC, RHS] = zeroDOF(CC, DD, RHS, bc, zeros(2,N));
        [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, zeros(2,N));
        RHS = RHS.';
        
        AA = AA(1:end-2,3:end);
        BB = BB(1:end-2,3:end);
        CC = CC(1:end-2,3:end);
        DD = DD(1:end-2,3:end);
        RHS = RHS(1:end-2,1:end-2);
        
        X22 = bartelsStewart(AA, BB, CC, DD, RHS );
        X12 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X22);
        X21 = ( bc(1:2,1:2).' \ (-bc(1:2,3:end)*X22.') ).';
        X11 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X21);
        X = [ X11 X12 ; X21 X22 ]; 
        Fapprox(:, :, k+(N+1)/2) = X;
    
    end

end



function [C1, E] = zeroDOF(C1, C2, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

    for ii = 1:size(B, 1) % For each boundary condition, zero a column.
        for kk = 1:size(C1, 1)
            if ( abs(C1(kk,ii)) > 10*eps )
                c = C1(kk, ii); % Constant required to zero entry out.
                C1(kk,:) = C1(kk,:) - c*B(ii,:);
                E(kk,:) = E(kk,:) - c*G(ii,:)*C2.';
            end
        end
    end

end


