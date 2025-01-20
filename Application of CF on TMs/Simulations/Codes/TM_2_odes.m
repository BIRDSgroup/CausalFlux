function dFdt = TM_2_odes(t,F)

%%%%  same as TM_3_New

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

R2 = ((F(10)+F(11)-(F(10)*F(11)))*F(1)*K)/(F(1)+Km);

R3 = (F(18)*F(1)*K)/(F(1)+Km);

R4 = (F(19)*F(2)*K)/(F(2)+Km);

R5 = ((F(19)+F(12)-(F(19)*F(12)))*F(3)*K)/(F(3)+Km);

R6 = ((F(14)+F(15)-(F(14)*F(15)))*F(3)*K)/(F(3)+Km);

R7 = ((F(21)+F(22)-(F(21)*F(22)))*K*F(4)*F(7))/(Km+F(4)+F(7)+(F(4)*F(7)));

R8 = ((F(15)+F(16)-(F(15)*F(16)))*F(5)*K)/(F(5)+Km);

R9 = (F(15)*F(5)*K)/(F(5)+Km);

R10 = Ks*F(6);

R11 = Ks*F(8);

% Metabolites in MN 
dFdt(1) = R1-R2-R3; % m1
dFdt(2) = R3-R4; % m4
dFdt(3) = R2+R4-R5-R6; % m2
dFdt(4) = R5-R7; % m3
dFdt(5) = R6-R8-R9; % m5 
dFdt(6) = R9-R10; % m6
dFdt(7) = R8-R7; % m7
dFdt(8) = R7-R11; % m8


%Genes in GRN
if co == "A" 
  dFdt(9) = 0;
else
  dFdt(9) = ti^(-1)*(Fact(F(8),1)-F(9)); % A
end 

dFdt(10) = ti^(-1)*(Fact(F(9),0)-F(10)); % B
dFdt(11) = ti^(-1)*(Fact(F(9),0)-F(11)); % C
dFdt(12) = ti^(-1)*(Fact(F(9),0)-F(12)); % D

if co == "E" 
  dFdt(13) = 0; % E
else
  dFdt(13) = ti^(-1)*(Const2-F(13)); % E
end 

dFdt(14) = ti^(-1)*(Fact(F(13),0)-F(14)); % F
dFdt(15) = ti^(-1)*(Fact(F(13),0)-F(15)); % G
dFdt(16) = ti^(-1)*(Fact(F(13),0)-F(16)); % H

if co == "E" 
  dFdt(17) = 0;
else
  dFdt(17) = ti^(-1)*(Const3-F(17)); % I
end 

dFdt(18) = ti^(-1)*(Fact(F(17),0)-F(18)); % J
dFdt(19) = ti^(-1)*(Fact(F(17),0)-F(19)); % K

if co == "E" 
  dFdt(20) = 0;
else
  dFdt(20) = ti^(-1)*(Fact(F(6),1)-F(20)); % X
end

dFdt(21) = ti^(-1)*(Fact(F(20),0)-F(21)); % L
dFdt(22) = ti^(-1)*(Fact(F(20),0)-F(22)); % M



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

