function [objfcn,electricityGenerationCost,electricityLoadCurtailmentCost,...
    gasPurchasingCost,gasCurtailmentCost] = objfcn_IEGSoperatingCost(Pg,LCe,PGs,Qptg,LCg,mpc,CDFe,CDFg,id,iGd)
Pg = mpc.baseMVA*Pg;
LCe = mpc.baseMVA*LCe;
%% 
electricityGenerationCost = sum( (Pg').^2 * mpc.gencost(:,5) + Pg' * mpc.gencost(:,6) + mpc.gencost(:,7) );
% electricityLoadCurtailmentCost = sum(LCe) * CDFe;
gasPurchasingCost = sum(PGs' * mpc.Gcost);
% gasCurtailmentCost = sum(LCg) * CDFg;
PTGsubsidy = sum(sum(Qptg))  *1e6/24 * 0.089 * 2.2/ 6.7;
PTGsubsidy = 0;
% 防止一直只切一个节点
gasCurtailmentCost = sum(LCg) * CDFg;
electricityLoadCurtailmentCost = sum(LCe) * CDFe;

objfcn = electricityGenerationCost + electricityLoadCurtailmentCost + gasPurchasingCost + gasCurtailmentCost - 1*PTGsubsidy;


end