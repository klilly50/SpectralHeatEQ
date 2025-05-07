classdef Operators
    properties
        C01 % Convert Cheb-T to Cheb-U
        C12 % Convert Cheb-U to C^2
        C02 % Convert Cheb-T to C^2
        D0 % Differentiate Cheb-T once
        D1 % Differentiate Cheb-T once and convert to Cheb-U
        D2 % Differentiate Cheb-T twice and convert to C^2
        R0 % Multiply Cheb-T by active variable
        R1 % Multiply Cheb-U by active variable
        R2 % Multiply C^2 by active variable
        A % R2^2*D2 + R2*C12*D1. Used in Laplacian.
        coF % R2^2*C02. Used in Laplacian.
        Dt % Differentiate Fourier once
        DCT % Values at Chebyshev nodes to Cheb-T coefficients
        DFT % Values at equidistant points to Fourier coefficients
        iDCT % Cheb-T to values at Chebyshev nodes
        iDFT % Fourier to values at equidistant points
        CC % Integrate Cheb-T over [-1, 1]
        CC2 % Integrate Cheb-T over [0, 1]
        CCf % Integrate Fourier over [0, 2*pi]
    end
    methods
        function obj = Operators(N,arg)
            if ( ~isnumeric(N) )
                error('Must be an integer value!')
            elseif ( floor(N) ~= N )
                error('Must be an integer value!')
            else
                if ( nargin < 2 )
                    arg = "C01, C12, C02, D0, D1, D2, R0, R1, R2, A, coF, Dt, DCT, DFT, iDCT, iDFT, CC, CC2, CCf";
                else
                    arg = string(arg);
                end
                
                if ( contains(arg,["C01" "C02" "D0" "coF" "D0"]) ) % C01
                    obj.C01 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); 
                    obj.C01(1, 1) = 1; 
                end
                if ( contains(arg,["C12" "C02" "R1" "D0" "A" "coF"]) ) % C12
                    K = 1./(1:N)';
                    obj.C12 = spdiags([K -K], [0 2], N, N);
                end
                if ( contains(string(arg),["C02" "coF"]) ) % C02
                    obj.C02 = obj.C01*obj.C12;
                end
                
                if ( contains(string(arg),["D1" "A" "D0"]) ) % D1
                    obj.D1 = spdiags((0:N)', 1, N, N);
                end
                if ( contains(string(arg),["D2" "A"]) ) % D2
                    obj.D2 = spdiags(2*(0:N)', 2, N, N);
                end
                if ( contains(string(arg),"D0") ) % D0
                    obj.D0 = obj.C01\obj.D1;
                end
                
                if ( contains(string(arg),"R0") ) % R0
                    obj.R0 = spdiags([1 .5; .5*ones(N,1) .5*ones(N,1)],[-1,1],N,N);
                end
                if ( contains(string(arg),["R2" "R1" "coF" "A"]) ) % R2
                    K = (1:N)'./(4:2:2*N+2)';
                    K1 = (3:N+2)'./(4:2:2*N+2)';
                    obj.R2 = spdiags([K K1], [-1 1], N, N);
                end
                if ( contains(string(arg),"R1") ) % R1
                    obj.R1 = obj.C12\obj.R2*obj.C12;
                end
                
                if ( contains(string(arg),"A") ) % A
                    obj.A = obj.R2^2*obj.D2 + obj.R2*obj.C12*obj.D1;
                end
                if ( contains(string(arg),"coF") ) % coF
                    obj.coF = obj.R2^2*obj.C02;
                end
                
                if ( contains(string(arg),"Dt") ) % Dt
                    if mod(N,2)
                        obj.Dt = spdiags(1i*(-N/2:N/2)', 0, N+1, N+1); 
                    else
                        obj.Dt = spdiags(1i*(-N/2:N/2)', 0, N+1, N+1); 
                    end
                end
                
                if ( contains(string(arg),"iDCT") ) % iDCT
                    obj.iDCT = chebtech2.coeffs2vals(eye(N)); 
                end
                if ( contains(string(arg),"iDFT") ) % iDFT
                    obj.iDFT = trigtech.coeffs2vals(eye(N)); 
                end
                if ( contains(string(arg),"DCT") ) % DCT
                    obj.DCT = chebtech2.vals2coeffs(eye(N)); 
                end
                if ( contains(string(arg),"DFT") ) % DFT
                    obj.DFT = trigtech.vals2coeffs(eye(N)); 
                end
                
                if ( contains(string(arg),["CC" "CC2"]) ) % CC
                    obj.CC = [2 kron(2./(1-(2:2:N).^2),[0 1])];
                    obj.CC = obj.CC(:,1:N);
                end
                if ( contains(string(arg),"CC2") ) % CC2
                    obj.CC2 = kron(1./((4:4:N+3)-2),[0 1 0 -1]); % Integrate on [0,1]
                    obj.CC2 = obj.CC2(1:N);
                    obj.CC2 = obj.CC2 + .5*obj.CC;
                end
                if ( contains(string(arg),"CCf") ) % CCf
                    obj.CCf = [zeros(1,(N/2)) 2*pi zeros(1,(N/2))];
                end
            end
        end
    end
end