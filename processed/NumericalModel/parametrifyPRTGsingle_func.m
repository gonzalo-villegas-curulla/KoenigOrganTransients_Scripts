function [PRTg] = parametrifyPRTGsingle_func(Ta_pass, Tb_pass, Tc_pass, Td_pass, Tomega_pass, SIG, ValveRampInit, Tend, tvec)

    PRTg = []; 
    if length(Ta_pass)>1
        for idx = 1 : length(Ta_pass)
            PRTg = [PRTg, run_simulation_func( Ta_pass(idx),Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tb_pass)>1
        for idx = 1 : length(Tb_pass)
            PRTg = [PRTg, run_simulation_func( Ta_pass,Tb_pass(idx),Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tc_pass)>1
        for idx = 1 : length(Tc_pass)
            PRTg = [PRTg, run_simulation_func( Ta_pass,Tb_pass,Tc_pass(idx),Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Td_pass)>1
        for idx = 1 : length(Td_pass)
            PRTg = [PRTg, run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass(idx),SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tomega_pass)>1
        for idx = 1 : length(Tomega_pass)
            PRTg = [PRTg, run_simulation_func( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass(idx), tvec)];
        end
    end
