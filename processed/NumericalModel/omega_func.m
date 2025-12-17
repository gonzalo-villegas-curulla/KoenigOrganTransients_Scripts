function OM = omega_func(t_ode, ValveRampInit, ValveRampEnd)
% == Initialisation function Omegat(t) ======================

    if     t_ode<=ValveRampInit
        OM = 0.0;
    elseif ValveRampEnd<t_ode
        OM = 1.0;        
    else
        OM = 0.5 - 0.5*cos(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit));
    end      
