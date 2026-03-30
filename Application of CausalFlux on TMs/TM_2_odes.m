
function dFdt = TM_1_odes(t,F)

Const2 = 0.9; % for gene E % 0.8
Const3 = 0.9; % for gene Y % 0.5
Const4 = 0.9; % for gene Z 
Km = F(14); % Michealis menten Kmax (5)
K =[F(15),F(15),F(15),F(15)];% enzyme constant in Vmax replacement; Vmax = [E]*K; [100,100,100,100]
Ks = F(16) ; % sink rate constant for R6, R7 and R8 -> 0.1



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
ti = F(17);


%Metabolites in MN
R1 = Const1*F(4)-F(1);
R2 = (((F(8)*F(9))*K(1)*F(1))/(Km+F(1)));
R3 = (((F(8)+F(7)-F(8)*F(7))*K(2)*F(1))/(Km+F(1)));
R4 = (((F(7)+F(9)-F(7)*F(9))*K(3)*F(2))/(Km+F(2)));
R5 = (((F(8)+F(9)-F(8)*F(9))*K(4)*F(3))/(Km+F(3)));
R6 = Ks*F(4);
R7 = Ks*F(5);
R8 = Ks*F(2);


dFdt(1) = R1-R2-R3; % d(m1)/dt = R1-R2-R3
dFdt(2) = R2-R4-R8; %d(m2)/dt = R2-R4-R8
dFdt(3) = R3-R5; %d(m3)/dt = R3-R5
dFdt(4) = R4-R6;  %d(m4)/dt = R4-R6
dFdt(5) =  R4+R5-R7; %d(m5)/dt = R4+R5-R7

%Genes in GRN
if co == "A" 
  dFdt(6) = 0; % A KO
else
  dFdt(6) = ti^(-1)*(Fact(F(2),1)+Fact(F(12),0)+Fact(F(13),0)-(Fact(F(2),1)*Fact(F(12),0)*Fact(F(13),0))-F(6)); % A
end 

if co == "B" 
  dFdt(7) = 0; %B KO
else
  dFdt(7) = ti^(-1)*(Fact(F(11),0)-F(7)); %B
end 

 dFdt(8) = ti^(-1)*((Fact(F(6),0)+Fact(F(7),0)-(Fact(F(6),0)*Fact(F(7),0)))-F(8)); % C
 %dFdt(8) = 0;

 dFdt(9) = ti^(-1)*((Fact(F(7),0)+Fact(F(10),0)-(Fact(F(7),0)*Fact(F(10),0)))-F(9)); % D
% dFdt(9) = 0; % D KO

if co == "E" 
   dFdt(10) = 0; % E KO
else
  dFdt(10) = ti^(-1)*(Const2-Fact(F(10),0)); % E
end

dFdt(12) = ti^(-1)*(Const3-Fact(F(12),0)); % Y

if co == "Z" 
   dFdt(13) = 0; % Z KO
else
  dFdt(13) = ti^(-1)*(Const4-Fact(F(13),0)); % Z
end

if co == "X" 
  dFdt(11)= 0; 
else
 dFdt(11)= ti^(-1)*((Fact(F(4),1)+Fact(F(5),1)-(Fact(F(4),1)*Fact(F(5),1)))-F(11)); % X
end

dFdt(14) = 0;
dFdt(15) = 0;
dFdt(16) = 0;
dFdt(17) = 0;

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




