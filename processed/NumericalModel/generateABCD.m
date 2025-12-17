function [Ta, Tb, Tc, Td] = generateABCD(sigma, Pf_targ, Pg_targ, PRTg, Vg, Vf)

    % From itemized recipe: 
    % * chose one value of $\Sigma$
    % * estimate A from PRTg
    % * deduced B,C,D from steady states
    % * simulate with very short $\tau_{\Omega}$ 


% Parse vars dimensions for loops and vectorized cases
sigma   = sigma(:);
Pf_targ = Pf_targ(:);
Pg_targ = Pg_targ(:);
PRTg    = PRTg(:);
Vg      = Vg(:);
Vf      = Vf(:);

% (1) -> Impose chosen sigma

% (2) 
Ta = PRTg/1.45; % 1.3, 1.4438, 1.45

% (3) Derive rest of characteristic times
    ABchunk = (1-Pg_targ)./(Pg_targ*(1-sigma^2) - Pf_targ + sigma^2);
    CDchunk = sqrt(Pf_targ ./ (Pg_targ*(1-sigma^2)-Pf_targ + sigma^2));
Tb = Ta./sqrt(ABchunk);
Tc = Tb.*Vf./Vg; 
Td = Tc.*CDchunk;

% (4) Omega_characteristic time: Def. 1ms
