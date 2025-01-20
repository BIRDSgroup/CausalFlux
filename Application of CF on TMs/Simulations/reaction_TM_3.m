function R = reaction_TM_3(F,ppK,ppV,ppkr,mulfac)

sf = mulfac;

Const1 = 3.2*sf; % for reaction R1 10.5 mmol/gDCW/hr

% Const2 = 0; % for gene E
Km = ppK; % Michealis menten Kmax 
K = ppV;  % enzyme constant in Vmax replacement; Vmax = [E]*K
Ks = ppkr ; % sink rate constant for R6, R7 and R8

R = zeros(size(F,1),9);
%Reactions

R(:,1) = Const1.*F(:,8)-F(:,1);


%R(:,2) = ((F(:,13)+F(:,14)+F(:,19)-(F(:,13).*F(:,14).*F(:,19))).*F(:,1).*K)./(F(:,1)+Km); % OR case

%%% new
R(:,2) = ((F(:,13)+F(:,14)-(F(:,13).*F(:,14))).*F(:,1).*K)./(F(:,1)+Km);

%R(:,3) = (F(:,15).*F(:,1).*K)./(F(:,1)+Km); 

%%% new
R(:,3) = ((F(:,15)+F(:,19)-(F(:,15).*F(:,19))).*F(:,1).*K)./(F(:,1)+Km);

R(:,4) = (F(:,16).*F(:,17).*F(:,2).*F(:,3).*K)./(F(:,2)+F(:,3)+(F(:,2).*F(:,3))+Km);

R(:,5) = (F(:,18).*F(:,3).*K)./(F(:,3)+Km);

%R(:,6) = (F(:,19).*F(:,4).*K)./(F(:,4)+Km);

%%% new 
R(:,6) = ((F(:,19)+F(:,22)-(F(:,19).*F(:,22))).*F(:,4).*K)./(F(:,4)+Km);

%R(:,7) = (F(:,20).*F(:,5).*F(:,4).*K)./(F(:,5)+F(:,4)+(F(:,5).*F(:,4))+Km);

%%% new
R(:,7) = ((F(:,20)+F(:,22)-(F(:,20).*F(:,22))).*F(:,5).*F(:,4).*K)./(F(:,5)+F(:,4)+(F(:,5).*F(:,4))+Km);

%R(:,8) = ((F(:,21)+F(:,22)-(F(:,21).*F(:,22))).*F(:,6).*F(:,7).*K)./(F(:,6)+F(:,7)+(F(:,6).*F(:,7))+Km); % OR case

%%% new
R(:,8) = ((F(:,21)+F(:,20)-(F(:,21).*F(:,20))).*F(:,6).*F(:,7).*K)./(F(:,6)+F(:,7)+(F(:,6).*F(:,7))+Km);

R(:,9) = Ks.*F(:,8);

end