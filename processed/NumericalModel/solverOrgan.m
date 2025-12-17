function [dydt] = solverOrgan(t_ode, y,A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

    omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);    
    dydt = zeros(2,1);
    dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A   -   real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
    dydt(2) =      real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C  -  real(sqrt(2*y(2)))/D;
