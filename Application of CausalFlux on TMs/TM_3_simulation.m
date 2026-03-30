curr_wd = 'D:\work\Integrated_network_model\TMS_actual_simulations\';
cd(curr_wd)
IC_DF = "IC_DF_TM_3_V1.csv";
IC_DF = readtable(IC_DF);

% Choose the Fic row index
Fic = IC_DF(1,:);


p_100 = "Params_100_km_vm_kr_t.csv";
p_100 = readtable(p_100);



para_100.K = p_100.Kms(1:100) ;para_100.V = p_100.Vms(1:100) ; para_100.kr = p_100.krs(1:100) ; para_100.t = p_100.t(1:100) ;


Cond = readvars("Condition.csv");

eb_var = readvars("Exch_b.csv");

co = Cond{2};

if eb_var == 3.2
    mulfac = 1;
elseif eb_var == 320
    mulfac = 100;
elseif eb_var == 3200
    mulfac = 1000;
end


tspan = 0:1:10000;
X = [0,0,0,0,0,0,0,0,0,0,0];
G = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
for i = (1:100)
  
    if co == "WT"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,Fic.F_14,Fic.F_15,Fic.F_16,Fic.F_17,Fic.F_18,Fic.F_19,Fic.F_20,Fic.F_21,Fic.F_22,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% WT
    elseif co == "X"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,Fic.F_14,Fic.F_15,Fic.F_16,Fic.F_17,Fic.F_18,Fic.F_19,0,Fic.F_21,Fic.F_22,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% X
    elseif co == "I"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,Fic.F_14,Fic.F_15,Fic.F_16,0,Fic.F_18,Fic.F_19,Fic.F_20,Fic.F_21,Fic.F_22,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% I
    elseif co == "A"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,0,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,Fic.F_14,Fic.F_15,Fic.F_16,Fic.F_17,Fic.F_18,Fic.F_19,Fic.F_20,Fic.F_21,Fic.F_22,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% A

    end

      [t,F_ako] = ode23s('TM_2_odes',tspan ,F0); % METHOD 1

      R_ = reaction_TM_2(F_ako,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,mulfac);
        
      G = [G;F_ako(end,9:22)];

      X = [X;R_(end,:)];
      
  
end


curr_wd = 'D:\work\Integrated_network_model\TMS_actual_simulations\';
cd(curr_wd)
%cd(curr_wd)
%writematrix(,'_A.csv')


if co == "WT"
      writematrix(G,['GE_simu_',eb_var,'.csv'])
end


writematrix(X,['ODE_simu_',co,'.csv'])




