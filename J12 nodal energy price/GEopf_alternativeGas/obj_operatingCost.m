function [f,totalCost,electricityGenerationCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg, mpc,CDF)
% unit: $/hour
Pg = mpc.baseMVA*Pg;
%

    electricityGenerationCost = sum( (Pg).^2 .* mpc.gencost(:,5) + Pg .* mpc.gencost(:,6) + mpc.gencost(:,7) );
%     electricityGenerationCost = sum(  Pg .* mpc.gencost(:,6)  + mpc.gencost(:,7) );
%     sum(totcost_yalmip(mpc.gencost, Pg(:,:)')); % Pg includes GFU, but the gencost are 0.
    gasPurchasingCost = sum(PGs' * mpc.Gcost);
    gasCurtailmentCost = sum(LCg) * CDF.gas;
% subsidy of hydrogen and methane productions
% totalCost = totalCost- 0.1/6*1e6*sum(sum(Qptg));  % original
PTGsubsidy = sum(sum(Qptg))  *1e6/24 * 0.089 * 2.2/ 6.7;  


%%
totalCost = electricityGenerationCost + gasPurchasingCost + gasCurtailmentCost - 0*PTGsubsidy;
f  = electricityGenerationCost + gasPurchasingCost + gasCurtailmentCost - 0*PTGsubsidy;
end