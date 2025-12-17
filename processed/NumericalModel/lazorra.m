clc;
clear;


MX1 = nan*ones(40,22);

for IIDX = 1 : 22
    for JDX = 1 : 40
        

        try
            sample_select = IIDX;
            TRANSNUM = JDX;

            run Model12_Transient_Meas_vs_Simul_UsesGenerateABCD_Fig09.m
            MX1(JDX,IIDX) = tauflow_simu;


        end

    end
end
