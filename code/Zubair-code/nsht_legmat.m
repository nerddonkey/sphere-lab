function [P Sc] = nsht_legmat(thetas, L, m)
%#codegen
% nsht_legmat - Compute Legendre matrix
%
% Compute the Legendre matrix.  Usage is given by
%
%  [P Sc] = nsht_mat(thetas, L, m)
%
% where thetas is the vector of theta samples, L is the harmonic band-limit
% and m is the order considered.  The computed matrix P is scaled version 
% and is ordered as P(ell, theta) where ell \in [m, L-1] and theta \in
% thetas. Sc is the scaling matrix such that P_ell(theta) = P(ell,
% theta).*10.^Sc.
%
% Notes:
%  - Kostelec recusrsion is implemented to compute the scaled legendre
%  coefficients for each theta \in thetas
%  - This code is converted to C for fast computation of Legendre matrices
%  and efficient implementation of forwad and inverse spherical harmonic
%  transform
%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details
%%
% Check arguments. 
if ( m >= L) 
  error('Require m < L.');
end

% following check is required if this function is only called by
% nsht_forward and nsht_inverse. We also call this function nsht_sharmonic
% and nsht_sharmonic_matrix.
% 
%if ( length(thetas) ~= L) 
%  error('Require L theta samples.');
%end

P = zeros(L-m, length(thetas));

Sc = zeros(size(P)); % scaling matrix contains exponents of 10


[Km power_10_first]  = norm_factor_first(m); % find scaling factor
    



for ii=1:1:length(thetas)
    TT = thetas(ii);
    if TT==0 && m~=0
            for ell=m:1:L-1
                P(ell-m+1,ii) = 0;
               Sc(ell-m+1,ii) =  1;
            end
    else

    [factor_remain2 power_10]  = sine_factor(TT, m);
    temp = power_10_first + power_10;
    dlm = Km*factor_remain2*(-1)^m*sqrt((2*m+1)/2)/sqrt(2*pi);

    dlm_1 = zeros(size(TT));

    for ell=m:1:L-1

        if abs(dlm)<10^(-2)
            while abs(dlm)<10^(-2)
                temp = temp-2;
                dlm = dlm*10^2;
                dlm_1 = dlm_1*10^2;
            end
        end
        if abs(dlm)>10^(2)
            while abs(dlm)>10^(2)
                temp = temp+2;
                dlm = dlm/10^2;
                dlm_1 = dlm_1/10^2;
            end
        end
        P(ell-m+1,ii) = dlm;

        Sc(ell-m+1,ii) =  temp;


        [ dlm1 ] = Kostelec_recursion_scaled(dlm, dlm_1,ell,m,TT );


        dlm_1 = dlm;
        dlm = dlm1;
    end

    end % end if TT==0
end



end


function [ dlm1 ] = Kostelec_recursion_scaled(dlm, dlm_1,ell,m,theta )

if (ell==0)
dlm1 = sqrt((2*ell+3)/(2*ell+1))* ((ell+1)*(2*ell+1)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(cos(theta).*dlm);
else
dlm1 = sqrt((2*ell+3)/(2*ell+1))* ((ell+1)*(2*ell+1)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(cos(theta).*dlm) - ...
       sqrt((2*ell+3)/(2*ell-1))* (sqrt((ell^2-m^2)*ell^2)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(ell+1)*dlm_1/ell; 
end
end


function [Km power_10_first]  = norm_factor_first(m)

Km=1;

power_10_first=0;


for ii=1:1:m
    
   Km = Km*sqrt((2*m-ii+1)/(m-ii+1))/2; 
    if abs(Km)>10^(2)
        while (abs(Km)>10^(2))
            power_10_first = power_10_first+2;
            Km = Km*10^(-2);
        end
    end
    if abs(Km)<10^(-2)
        while abs(Km)<10^(-2)
            power_10_first = power_10_first-2;
            Km = Km*10^(2);
        end
    end
    
end


end


function [factor_remain power_10]  = sine_factor(TT, m)
if TT==0 && m~=0
   factor_remain=0;
   power_10=0;
else
factor_remain = 1;
power_10=0;
    for ii=1:1:m
        factor_remain = factor_remain*sin(TT);
        if abs(factor_remain)<10^(-2)
            while abs(factor_remain)<10^(-2)
                power_10 = power_10-2;
                factor_remain = factor_remain*10^2;
            end
        end
    end
end
end