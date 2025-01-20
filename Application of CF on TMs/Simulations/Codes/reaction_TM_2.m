function R = reaction_TM_2(F,ppK,ppV,ppkr, mulfac)

sf = mulfac;

Const1 = 3.2*sf; % for reaction R1 10.5 mmol/gDCW/hr

% Const2 = 0; % for gene E
Km = ppK; % Michealis menten Kmax 
K = ppV;  % enzyme constant in Vmax replacement; Vmax = [E]*K
Ks = ppkr ; % sink rate constant for R6, R7 and R8

R = zeros(size(F,1),11);
%Reactions
R(:,1) = Const1.*F(:,8)-F(:,1);

R(:,2) = ((F(:,10)+F(:,11)-(F(:,10).*F(:,11))).*F(:,1).*K)./(F(:,1)+Km);

R(:,3) = (F(:,18).*F(:,1).*K)./(F(:,1)+Km);

R(:,4) = (F(:,19).*F(:,2).*K)./(F(:,2)+Km);

R(:,5) = ((F(:,19)+F(:,12)-(F(:,19).*F(:,12))).*F(:,3).*K)./(F(:,3)+Km);

R(:,6) = ((F(:,14)+F(:,15)-(F(:,14).*F(:,15))).*F(:,3).*K)./(F(:,3)+Km);

R(:,7) = ((F(:,21)+F(:,22)-(F(:,21).*F(:,22))).*F(:,4).*F(:,7).*K)./(Km+F(:,4)+F(:,7)+(F(:,4).*F(:,7)));

R(:,8) = ((F(:,15)+F(:,16)-F(:,15).*F(:,16)).*F(:,5).*K)./(F(:,5)+Km);

R(:,9) = (F(:,15).*F(:,5).*K)./(F(:,5)+Km);

R(:,10) = Ks.*F(:,6);

R(:,11) = Ks.*F(:,8);

end