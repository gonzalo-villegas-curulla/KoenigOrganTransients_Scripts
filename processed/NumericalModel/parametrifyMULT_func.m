function [PRTg, PRTf] = parametrifyMULT_func(Ta_pass,...
    Tb_pass, ...
    Tc_pass, ...
    Td_pass, ...
    Tomega_pass, ...
    SIG, ...
    ValveRampInit, ...
    Tend, ...
    tvec)

    PRTg = [];
    PRTf = []; 
    if length(Ta_pass)>1
        for idx = 1 : length(Ta_pass)
            [tmp(1),tmp(2)] = run_simulation_func( Ta_pass(idx),Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);            
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if length(Tb_pass)>1
        for idx = 1 : length(Tb_pass)
            [tmp(1),tmp(2)] = run_simulation_func( Ta_pass,Tb_pass(idx),Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if length(Tc_pass)>1
        for idx = 1 : length(Tc_pass)
            [tmp(1),tmp(2)] = run_simulation_func( Ta_pass,Tb_pass,Tc_pass(idx),Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if length(Td_pass)>1
        for idx = 1 : length(Td_pass)
            [tmp(1),tmp(2)] = run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass(idx),SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if length(Tomega_pass)>1
        for idx = 1 : length(Tomega_pass)
            [tmp(1),tmp(2)] = run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass(idx), tvec);
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if length(SIG)>1
        for idx = 1 : length(SIG)
            [tmp(1), tmp(2)] = run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG(idx),     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);
            PRTg = [PRTg, tmp(1)];
            PRTf = [PRTf, tmp(2)];
        end
    end
    if ( (length(Ta_pass)==1)  && (length(Tb_pass)==1)  && (length(Tc_pass)==1)  && (length(Td_pass)==1)  && (length(Tomega_pass)==1)  && (length(SIG)==1) )
        [PRTg, PRTf] = run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec);        
        fprintf('Using last opt of parametrify func\n');
    end

