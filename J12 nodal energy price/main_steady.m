clear
clc
yalmip('clear')
load("windData.mat")
%% parameters
[mpc,gtd] = case24GE_J15(); 

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
baseMVA = 100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);

nb   = size(mpc.bus, 1);    %% number of buses
ng = size(mpc.gen,1);
nd = size(find(mpc.bus(:,PD)~=0),1);
id = find(mpc.bus(:,PD)~=0);
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl = size(mpc.Gline,1);
nGen = sum(mpc.gen(:,22)==1)+sum(mpc.gen(:,22)==0);% all gen(TFU and GFU), excluded dispatchable loads
nGPP = size(mpc.gfuIndex,1);
nGs = size(mpc.Gsou,1);
nGd = size(find(mpc.Gbus(:,3)~=0),1);
nPTG = size(mpc.ptg,1);
iGd = find(mpc.Gbus(:,3)~=0);
nWindFarm = size(mpc.windfarmIndex,1);

mpc0 = mpc;
mpc.gen(:,PMIN) = 0; % 最小出力为0
mpc.Gbus(:,3) = mpc.Gbus(:,3)*46.0980/48.9660; % 减少负荷数量，不然总气源能量不足
% mpc.branch([28,31,32,33],6) = mpc.branch([28,31,32,33],6)/2; % 围绕着RNG的线路容量/2
mpc.Gbus(:,6) = 99;%气压上限全部增加
mpc.Gbus(:,5) = 25;%气压下限全部变一致
mpc.Gsou(:,3) = 0;

    % ------test ------
%     mpc.gen(9:11,PMAX) = 80;
%     mpc.gencost(9:11,5:7) = 20 * mpc.gencost(9:11,5:7);
%     mpc.GEcon=[
%     4 15
%     14 1
%     16 2
%     1 7
% ];
    %------------
% GEopf_alternativeGas_nonlinear(mpc);
% [information_milp, GEresult_milp] = GEopf_milp(mpc);
% ------case1-------------
% mpc.branch(:,6) = mpc.branch(:,6)/1.3; 
% % tic
% [information, solution] = GEopf_alternativeGas(mpc,[]);
% toc
%--------------------------



NK = 48;

% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
nGasType = 7; iCombustibleGas = 1:5; iNonCombustibleGas = 6:7;
[GCV, M, fs, a, R, T_stp, Prs_stp, Z_ref, T_gas, eta, CDF] = initializeParameters_J15();

%% wind data and steady-state results
% windSpeed0 = windData(:,1:nWindFarm); %如果用平均值，会损失波动性，不如随便挑一天
% for iWindFarm = 1:nWindFarm
%     windSpeed(:,iWindFarm) = interp(windSpeed0(:,iWindFarm),2);
% end
% for k = 1:NK
%     for iWindFarm = 1:nWindFarm
%         windTurbinePower(k,iWindFarm) = windTurbine(windSpeed(k,iWindFarm),mpc.gen(mpc.windfarmIndex(iWindFarm),PMAX));
%     end
% end
% [gasDemandCurve0,ElectricityDemandCurve0] = energyDemandProfile();
% gasDemandCurve = interp(gasDemandCurve0,2); ElectricityDemandCurve = interp(ElectricityDemandCurve0,2);
% kg = gasDemandCurve/max(gasDemandCurve); ke = ElectricityDemandCurve/max(ElectricityDemandCurve);
% ------ case 2 ------:
% 设定风机出力的边界条件
% mpc.branch(:,6) = mpc.branch(:,6)/2; 
% for k = 1:NK
%     for iWindFarm = 1:nWindFarm
%         windTurbinePower(k,iWindFarm) = windTurbine(windSpeed(k,iWindFarm),mpc.gen(mpc.windfarmIndex(iWindFarm),PMAX));
%     end
% end
% for k = 1:NK
%     yalmip('clear')
%     newmpc = mpc;
% 
%     newmpc.gen(mpc.windfarmIndex,PMAX) = windTurbinePower(k,:)';
%     newmpc.Gbus(:,3) = mpc.Gbus(:,3) * kg(k); 
%     newmpc.bus(:,PD) = mpc.bus(:,PD) * ke(k);
%     [information{k}, solution{k}] = GEopf_alternativeGas(newmpc,[]);
%     % results
%     nodalElectricityPrice(k,:) = solution{k}.nodalElectricityPrice';
%     nodalGasPrice(k,:) = solution{k}.nodalGasPrice';
% 
% end

%% case 2.2
% mpc.branch(:,6) = mpc.branch(:,6)/2; 
% for k = 1:NK
%     for iWindFarm = 1:nWindFarm
%         windTurbinePower(k,iWindFarm) = windTurbine(windSpeed(k,iWindFarm),mpc.gen(mpc.windfarmIndex(iWindFarm),PMAX));
%     end
% end
% for k = 1:NK
%     yalmip('clear')
%     newmpc = mpc;
% 
%     newmpc.gen(mpc.windfarmIndex,PMAX) = windTurbinePower(k,:)';
%     newmpc.Gbus(:,3) = mpc.Gbus(:,3) * kg(k); 
%     newmpc.bus(:,PD) = mpc.bus(:,PD) * ke(k);
% 
%     options.FSlimit = 0.2;
%     [information1{k}, solution1{k}] = GEopf_alternativeGas(newmpc,options);
%     nodalElectricityPrice1(k,:) = solution1{k}.nodalElectricityPrice';
%     nodalGasPrice1(k,:) = solution1{k}.nodalGasPrice';
%     PTGproduction1(k,:) = sum(solution1{k}.Qptg,1);
%     hyCompsition1(k,:) = solution1{k}.gasComposition(:,5)';
% 
%     options.FSlimit = 0.1;
%     [information2{k}, solution2{k}] = GEopf_alternativeGas(newmpc,options);
%     nodalElectricityPrice2(k,:) = solution2{k}.nodalElectricityPrice';
%     nodalGasPrice2(k,:) = solution2{k}.nodalGasPrice';
%     PTGproduction2(k,:) = sum(solution2{k}.Qptg,1);
%     hyCompsition2(k,:) = solution2{k}.gasComposition(:,5)';
% end
% impact.gasComposition = (hyCompsition1 - hyCompsition2)./(hyCompsition1+hyCompsition2);
% impact.nodalElectricityPrice = (nodalElectricityPrice1 - nodalElectricityPrice2)./(nodalElectricityPrice1+nodalElectricityPrice2);
% impact.nodalGasPrice = (nodalGasPrice1 - nodalGasPrice2)./(nodalGasPrice1+nodalGasPrice2);

%% case 3
% load("electricityGasLoadData_Australia.mat")
% mpc.branch(:,6) = mpc.branch(:,6)/2; 
% ke3 = electricityLoad./max(electricityLoad); 
% kg3 = gasLoad./max(gasLoad);
% K_notConverge = [];
% for k = 1:8760*2
% % for k = 2913
%     % wind
%     for iWindFarm = 1:nWindFarm
%             windTurbinePower(k,iWindFarm) = windTurbine(windSpeed(k,iWindFarm),mpc.gen(mpc.windfarmIndex(iWindFarm),PMAX));
%     end
%         yalmip('clear')
%         newmpc = mpc;
%     
%         newmpc.gen(mpc.windfarmIndex,PMAX) = windTurbinePower(k,:)';
%         newmpc.Gbus(:,3) = mpc.Gbus(:,3) * kg3(k); 
%         newmpc.bus(:,PD) = mpc.bus(:,PD) * ke3(k);
%         [information{k}, solution{k},solverTime(k)] = f_GEopf_alternativeGas(newmpc,[]);
%         if information{k}.ruleoutFlag == 1 % not converge
%             K_notConverge = [K_notConverge; k];
%         end
%         % results
%         nodalElectricityPrice(k,:) = solution{k}.nodalElectricityPrice';
%         nodalGasPrice(k,:) = solution{k}.nodalGasPrice';
%     
% end
% k=2913
%% results
% nodalElectricityPrice1 = zeros(NK,nb,365); nodalGasPrice1 = zeros(NK,nGb,365);
% for k = 1:8760*2
%     iDay = floor((k-1)/NK) + 1;
%     ik = k-(iDay-1)*NK;
%     ke31(ik,iDay) = ke3(k);
%     nodalElectricityPrice1(ik,:,iDay) = solution{k}.nodalElectricityPrice';
%     nodalGasPrice1(ik,:,iDay) = solution{k}.nodalGasPrice';
% end
% % 注意初始时间是8：00开始的
% % risk
% expectedNodalElectricityPrice = mean(nodalElectricityPrice1,3);
% expectedNodalGasPrice = mean(nodalGasPrice1,3);
% stdNodalElectricityPrice = std(nodalElectricityPrice1,[],3);
% stdNodalGasPrice = std(nodalGasPrice1,[],3);
% 
% corElectricity = zeros(nb,nb); corGas = zeros(nGb,nGb);
% corElectricityAndGas = zeros(nb,nGb);
% for iDay = 1:365
%     % electricity
% for iEB = 1:nb
%     for jEB = 1:nb
%         [corElectricity(iEB,jEB,iDay)] = xcorr(nodalElectricityPrice1(:,iEB,iDay),nodalElectricityPrice1(:,jEB,iDay),0,'normalized');
%     end
% end
% %gas
% for iGB = 1:nGb
%     for jGB = 1:nGb
%         [corGas(iGB,jGB,iDay)] = xcorr(nodalGasPrice1(:,iGB,iDay),nodalGasPrice1(:,jGB,iDay),0,'normalized');
%     end
% end
% % electricity and gas
% for iEB = 1:nb
%     for jGB = 1:nGb
%         [corElectricityAndGas(iEB,jGB,iDay)] = xcorr(nodalElectricityPrice1(:,iEB,iDay),nodalGasPrice1(:,jGB,iDay),0,'normalized');
%     end
% end
% end
% corElectricity = mean(corElectricity,3);
% corGas = mean(corGas,3);
% corElectricityAndGas = mean(corElectricityAndGas,3);

%% sensitivity analysis
% mpc.branch(:,6) = mpc.branch(:,6)/2; 
% mpc.ptg(:,6) = 4;
% for i1 = 1:3
%     for i2 = 1:3
%         yalmip('clear')
%         
%         newmpc = mpc;
%         options.FSlimit = 0.1 + (i1-1)*0.1;
%         options.WIlimit = i2;
%     
%         [information{i1,i2}, solution{i1,i2}] = GEopf_alternativeGas(newmpc,options);
%         nodalElectricityPrice{i1,i2} = solution{i1,i2}.nodalElectricityPrice;
%         nodalGasPrice{i1,i2} = solution{i1,i2}.nodalGasPrice;
%         PTGproduction{i1,i2} = sum(solution{i1,i2}.Qptg,1);
%         hyCompsition{i1,i2} = solution{i1,i2}.gasComposition(:,5);
% 
%         nodalElectricityPriceAll(:,3*(i1-1)+i2) = nodalElectricityPrice{i1,i2};
%         nodalGasPriceAll(:,3*(i1-1)+i2) = nodalGasPrice{i1,i2};
%     end
% 
% end
% impactLowFSupWI = (nodalGasPrice{1,3} - nodalGasPrice{1,2}) ./ nodalGasPrice{1,2};
% impactLowFSdownWI = (nodalGasPrice{1,1} - nodalGasPrice{1,2}) ./ nodalGasPrice{1,2};
% impactMidFSupWI = (nodalGasPrice{2,3} - nodalGasPrice{2,2}) ./ nodalGasPrice{2,2};
% impactMidFSdownWI = (nodalGasPrice{2,1} - nodalGasPrice{2,2}) ./ nodalGasPrice{2,2};
% impactHighFSupWI = (nodalGasPrice{3,3} - nodalGasPrice{3,2}) ./ nodalGasPrice{3,2};
% impactHighFSdownWI = (nodalGasPrice{3,1} - nodalGasPrice{3,2}) ./ nodalGasPrice{3,2};
% 
% impactLowWIupFS = (nodalGasPrice{3,1} - nodalGasPrice{2,1}) ./ nodalGasPrice{2,1};
% impactLowWIdownFS = (nodalGasPrice{1,1} - nodalGasPrice{2,1}) ./ nodalGasPrice{2,1};
% impactMidWIupFS = (nodalGasPrice{3,2} - nodalGasPrice{2,2}) ./ nodalGasPrice{2,2};
% impactMidWIdownFS = (nodalGasPrice{1,2} - nodalGasPrice{2,2}) ./ nodalGasPrice{2,2};
% impactHighWIupFS = (nodalGasPrice{3,3} - nodalGasPrice{2,3}) ./ nodalGasPrice{2,3};
% impactHighWIdownFS = (nodalGasPrice{1,3} - nodalGasPrice{2,3}) ./ nodalGasPrice{2,3};
% 
% impactAnalysis = [impactLowFSupWI, impactLowFSdownWI, impactMidFSupWI, impactMidFSdownWI, impactHighFSupWI, impactHighFSdownWI, ...
%     impactLowWIupFS, impactLowWIdownFS, impactMidWIupFS, impactMidWIdownFS, impactHighWIupFS, impactHighWIdownFS];

%% large case
[mpc] = NorthWestCase();
% does not consider gas pressure drop
mpc.Gsou(:,3) = 0;
mpc.ptg(:,6) = 3;
% 把ptg往风机在的节点上调
% mpc.ptg(:,2) = [
%     111
%     106
%     1
%     124
%     147
%     163
%     174
%     196
% ];
% mpc.branch(:,6) = mpc.branch(:,6)/1.5;
% [information_milp, GEresult_milp] = GEopf_milp_largeCase(mpc);

[information, solution] = GEopf_alternativeGas_largeCase(mpc,[]);

save 