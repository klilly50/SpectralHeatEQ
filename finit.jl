# Import Packages!

using Pkg
Pkg.add("SpecialFunctions")

using Random
using SpecialFunctions

# Exact Solution Function u(r, z, th, t)

w54 = 18.98013387517992112;  # 4th root of the J-Bessel function of order 5
w33 = 13.0152007216984344;   # 3rd root of the J-Bessel function of order 5
q1 = 3;                      # Oscillation in z.
q2 = 2;                      

function u(r, z, th, t)

    u(r, z, th, t) = 10*exp(-((pi*q1)^2 + w33^2)*t)*sin(q1*pi*z).*besselj(3,w33*r).*sin(3*th) + exp(-((pi*q2)^2 + w54^2)*t)*sin(q2*pi*z).*besselj(5,w54*r).*cos(5*th);

    # alt option? u(r,z,th,t) = 10*exp(-((pi*q1)^2 + w33^2)*t)*sin(q1*pi*z).*besselj(3,w33*r).*sin(3*th);

# besselj computes bessel function of the first kind. What does this mean??
# julia equivalent is the same call out with special functions package

return

display(u(r, z, th, t))

end

# Evaluate at a n-by-n-by-n grid at time t=T*dt

dt = 1e-9; 
T = 200;
T_FD = 6;

finit(r,z,th) = u(r,z,th,0);
