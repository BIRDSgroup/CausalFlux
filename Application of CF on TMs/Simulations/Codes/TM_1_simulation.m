
curr_wd = 'D:\work\Integrated_network_model\TMS_actual_simulations\';
cd(curr_wd)
IC_DF = "IC_DF.csv";
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
X = [0,0,0,0,0,0,0,0];
G = [0,0,0,0,0,0,0,0];
for i = (1:100)
  
    if co == "WT"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% WT
    elseif co == "B"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,0,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)]; % B KO
    elseif co == "E"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,0,Fic.F_11,Fic.F_12,Fic.F_13,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)]; % E KO
    elseif co == "Z"
        F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,0,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)]; % Z KO
    elseif co == "A"
         F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,0,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,Fic.F_11,Fic.F_12,Fic.F_13,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)];% A KO
    elseif co == "X"
         F0 = [Fic.F_1,Fic.F_2,Fic.F_3,Fic.F_4,Fic.F_5,Fic.F_6,Fic.F_7,Fic.F_8,Fic.F_9,Fic.F_10,0,Fic.F_12,Fic.F_13,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,p_100.t(i)]; % X KO
    end

      [t,F_ako] = ode23s('TM_1_odes',tspan ,F0); % METHOD 1

      R_ = reaction_TM_1(F_ako,p_100.Kms(i),p_100.Vms(i)*mulfac,p_100.krs(i)*mulfac,mulfac);
        
      G = [G;F_ako(end,6:13)];

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




