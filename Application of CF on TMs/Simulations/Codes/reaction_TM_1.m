function R = reaction_TM_1(F,ppK,ppV,ppkr,mulfac)

sf = mulfac;

Const1 = 3.2*sf; % for reaction R1 10.5 mmol/gDCW/hr
% Const2 = 0; % for gene E
Km = ppK; % Michealis menten Kmax 
K =  [ppV,ppV,ppV,ppV];  % enzyme constant in Vmax replacement; Vmax = [E]*K
Ks =ppkr ; % sink rate constant for R6, R7 and R8

R = zeros(size(F,1),8);
%Metabolites in MN
R(:,1) = Const1.*F(:,4)-F(:,1);
R(:,2) = (((F(:,8).*F(:,9)).*K(1).*F(:,1))./(Km+F(:,1)));
R(:,3) = (((F(:,8)+F(:,7)-F(:,8).*F(:,7)).*K(2).*F(:,1))./(Km+F(:,1)));
R(:,4) = (((F(:,7)+F(:,9)-F(:,7).*F(:,9)).*K(3).*F(:,2))./(Km+F(:,2)));
R(:,5) = (((F(:,8)+F(:,9)-F(:,8).*F(:,9)).*K(4).*F(:,3))./(Km+F(:,3)));
R(:,6) = Ks*F(:,4);
R(:,7) = Ks*F(:,5);
R(:,8) = Ks*F(:,2);
end