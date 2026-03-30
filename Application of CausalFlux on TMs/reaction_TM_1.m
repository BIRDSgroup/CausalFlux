function dFdt = TM_3_odes(t,F)

Km = F(23); % Michealis menten Kmax (5)
K = F(24);% enzyme constant in Vmax replacement; Vmax = [E]*K; [100,100,100,100]
Ks = F(25) ; % sink rate constant for R6, R7 and R8 -> 0.1

curr_wd = 'D:\work\Integrated_network_model\TMS_actual_simulations\';
cd(curr_wd)
Cond = readvars("Condition.csv");
eb_var = readvars("Exch_b.csv");

if eb_var == 3.2
    mulfac = 1;
elseif eb_var == 320
    mulfac = 100;
elseif eb_var == 3200
    mulfac = 1000;
end

co = Cond{2};

Const1 = 3.2*mulfac; % for reaction R1 10.5 mmol/gDCW/hr


%ti = 100; %[100]
ti = F(26);
Const2 = 0.9;
Const3 = 0.9;

% Reactions (rates) in MN
R1 = Const1*F(8)-F(1);
 
R2 = ((F(13)+F(14)-(F(13)*F(14)))*F(1)*K)/(F(1)+Km);

R3 = ((F(15)+F(19)-(F(15)*F(19)))*F(1)*K)/(F(1)+Km);


R4 = (F(16)*F(17)*F(2)*F(3)*K)/(F(2)+F(3)+(F(2)*F(3))+Km);

R5 = (F(18)*F(3)*K)/(F(3)+Km);

R6 = ((F(19)+F(22)-(F(19)*F(22)))*F(4)*K)/(F(4)+Km);

R7 = ((F(20)+F(22)-(F(20)*F(22)))*F(5)*F(4)*K)/(F(5)+F(4)+(F(5)*F(4))+Km);

R8 = ((F(21)+F(20)-(F(21)*F(20)))*F(6)*F(7)*K)/(F(6)+F(7)+(F(6)*F(7))+Km);

R9 = Ks*F(8);

% Metabolites in MN 
dFdt(1) = R1-R2; %m1
dFdt(2) = R2-R4; %m2
dFdt(3) = R3-R4-R5; %m3
dFdt(4) = R4-R6-R7; %m4
dFdt(5) = R5-R7; %m5
dFdt(6) = R6-R8; %m6
dFdt(7) = R7-R8; %m7
dFdt(8) = R8-R9; %m8

%Genes in GRN
if co == "A" 
  dFdt(9) = 0;
else
  dFdt(9) = ti^(-1)*(Fact(F(8),1)-F(9)); % A
end 
 

dFdt(10) = ti^(-1)*(Const2-F(10)); % E
% dFdt(10) = 0;

  dFdt(11) = ti^(-1)*(Const3-F(11)); % I
%  dFdt(11) = 0;

if co == "X" 
  dFdt(12) = 0;
else
  dFdt(12) = ti^(-1)*(Fact(F(8),1)-F(12)); % X
end 
 
dFdt(13) = ti^(-1)*(Fact(F(9),0)-F(13)); % B
dFdt(14) = ti^(-1)*((Fact(F(9),0)+Fact(F(10),0)-(Fact(F(9),0)*Fact(F(10),0)))-F(14)); % C
dFdt(15) = ti^(-1)*((Fact(F(9),0)+Fact(F(10),0)-(Fact(F(9),0)*Fact(F(10),0)))-F(15)); % D

dFdt(16) = ti^(-1)*(Fact(F(10),0)-F(16)); % F
dFdt(17) = ti^(-1)*((Fact(F(11),0)+Fact(F(10),0)-(Fact(F(11),0)*Fact(F(10),0)))-F(17)); % G

dFdt(18) = ti^(-1)*(Fact(F(11),0)-F(18)); % H
dFdt(19) = ti^(-1)*((Fact(F(11),0)+Fact(F(12),0)-(Fact(F(11),0)*Fact(F(12),0)))-F(19)); % J
dFdt(20) = ti^(-1)*((Fact(F(11),0)+Fact(F(12),0)-(Fact(F(11),0)*Fact(F(12),0)))-F(20)); % K

dFdt(21) = ti^(-1)*(Fact(F(12),0)-F(21)); % L
dFdt(22) = ti^(-1)*(Fact(F(12),0)-F(22)); % M

dFdt(23) = 0;
dFdt(24) = 0;
dFdt(25) = 0;
dFdt(26) = 0;
% dFdt(27) = 0;
% dFdt(28) = 0;

dFdt = dFdt'; 
end

function F = Fact(C,bool)
if bool==0
    EC50 = 0.5; % fractional activation of input species required to induce half maximal activation of output species
    n = 1.4; % Hill coefficient
    beta = ((EC50^n)-1)/((2*EC50^n)-1);
    K = (beta-1).^(1/n);
    F = (beta*C^n)/(K^n + C^n);
elseif bool==1
    if C>=0.5
        F=1;
    else
        F=0;
    end
end
end






%%% Keeping a track of good params

% 0.1,50,5
