function [objfcn,constraints,GEopf_alternativeGas_optimizer] = GEopf_alternativeGas_construct(gasComposition_ref, Qd_ref, Qgpp_ref, gasFlow_sum_ref, newC, W_ref, alpha_PHI, alpha_x, alpha_Qd, ...
    mpc,nb,nGb,nGl,nGs,nGd,nGasType,nPTG,nGPP,iGd,iGppGd,signGf,REF,GCV,M,S_ng,FSlimit,FS_ng,GCV_ng,eta, ...
    fs,CDF,solvedFlag)
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
% 1. state variables ======================================================
Prs_square = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGd,nGasType);% Mm3/day
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = sdpvar(nPTG,1); % electricity consumption, 1/100 MW
Pg = sdpvar(size(mpc.gen,1),1); % include TPP, GPP and renewable generators, 1/100 MW
Qgpp = sdpvar(nGPP, nGasType); % Mm3/day
Va = sdpvar(nb,1);
gasComposition = sdpvar(nGb,nGasType); 
gasFlow = sdpvar(nGl,nGasType);% Mm3/day
sigma_PHI = sdpvar(nGl,1); % error limit for gas flow
sigma_x = sdpvar(nGb, nGasType); % error for gas composition
sigma_Qd = sdpvar(nGd, nGasType);
sigma_Qgpp = sdpvar(nGPP,nGasType);
sigma_Gf = sdpvar(nGl,nGasType);
% test
LCg = zeros(nGd,1);
gamma = signGf;
% 2. bounds ===============================================================
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6); % bar
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4); % Mm3/day
Qptgmin = mpc.ptg(:,5); Qptgmax = mpc.ptg(:,6);
QptgMax_hydrogen = Qptgmax;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
LCgmin = zeros(nGd,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0;  
refs = find(mpc.bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = 1;   %% voltage angle reference constraints
Val(refs) = 1;
gasFlowMax = mpc.Gline(:,5);
% 3. constraints ==========================================================
% gas source cons
gasSourceCons = [  PGsmin <= PGs <= PGsmax;];
% gas demand cons
energyDemand = mpc.Gbus(iGd,3) * GCV_ng; % energy need of these gas bus
gasDemandCons = [
     -sigma_Qd <= gasComposition_ref(iGd,:) .* repmat(sum(Qd_ref,2),[1,nGasType]) + (gasComposition(iGd,:)-gasComposition_ref(iGd,:)) .* repmat(sum(Qd_ref,2),[1,nGasType]) ...
        + gasComposition_ref(iGd,:) .* (repmat(sum(Qd,2),[1,nGasType])-repmat(sum(Qd_ref,2),[1,nGasType])) - Qd <= sigma_Qd; 
    Qd >= 0;
    sigma_Qd >= 0;
    ]:'gasDemandCons';
gasDemandBalanceCons = [    Qd * GCV.all'/1e9 == energyDemand/1e9;   ]:'gasDemandBalanceCons';
% nodal gas flow balance cons
nodalGasFlowBalanceCons = [
    consfcn_nodalGasFlowBalance(PGs,Qd,Qgpp,Qptg, gasFlow,mpc,nGasType,nGPP,nGd,iGd) == 0;
    ]:'nodalGasFlowBalanceCons';
% ptg
PTGcons = [
    ( Qptg(:,1) * 1e6/24/3600 *GCV.CH4 / eta.methanation + Qptg(:,2) * 1e6/24/3600 * GCV.hy ...
        ) /1e6 == Pptg*baseMVA / eta.electrolysis; % w
    0 <= Pptg*baseMVA / eta.electrolysis <= QptgMax_hydrogen /24/3600 * GCV.hy; % 如果全用来制氢，可制备xx
    0 <= Qptg;
    ];
% gpp
Pgpp = Pg(mpc.gfuIndex) * baseMVA;% MW
GPPcons = [
    -sigma_Qgpp <= gasComposition_ref(iGppGd,:) .* repmat(sum(Qgpp_ref,2),[1,nGasType]) + (gasComposition(iGppGd,:)-gasComposition_ref(iGppGd,:)) .* repmat(sum(Qgpp_ref,2),[1,nGasType]) ...
        + gasComposition_ref(iGppGd,:) .* (repmat(sum(Qgpp,2),[1,nGasType])-repmat(sum(Qgpp_ref,2),[1,nGasType])) - Qgpp <= sigma_Qgpp;
    sigma_Qgpp >= 0;
    Pgpp == Qgpp/24/3600 * GCV.all'*eta.GFU;
    Qgpp >= 0;
    ]:'GPPcons';
% electricity flow
electricityCons = [
    consfcn_electricBranchFlow(Va, mpc, il) <= 0;
    Pgmin <= Pg <= Pgmax;
    Va(refs) == 0; % 除了slackbus外，其他相角都没约束。但是一些solver在处理inf的上下限的时候有问题
    ]:'electricityCons';
electricityBalanceCons = [consfcn_electricPowerBalance(Va,Pg,Pptg,mpc) == 0;]:'electricityBalanceCons';
% wobbe index
GCV_nodal = gasComposition * GCV.all';
S_nodal = gasComposition * M.all' / M.air;
sqrtS = 0.5 * (S_nodal/sqrt(S_ng) + sqrt(S_ng));
WI_nodal = GCV_nodal ./ sqrtS; % 如果不能自动转化，那就手动化一下
WImin = 46.5; WImax = 52.85;
WobbeIndexCons = [
    WImin * sqrtS <= GCV_nodal/1e6 <= WImax * sqrtS
    ]:'WobbeIndexCons';
% Weaver flame speed factor
FSnodal = gasComposition * fs.All';
FScons = [
    (1-FSlimit) * FS_ng <= FSnodal <= (1+FSlimit) * FS_ng;
    ]:'FScons';
% specific gravity
SGcons = [
    S_nodal <= 0.7;
    ]:'SGcons';
% other security cons
otherSecurityCons = [
    0 <= gasComposition <= 1;
    ]:'otherSecurityCons';
% SOC reformulation for gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
gasFlow_sum = sum(gasFlow,2);
% PrsSquareLeft = (1-gamma)/2 .* Prs(FB).^2 + (1+gamma)/2 .* Prs(TB).^2; 
% PrsSquareRight = (1+gamma)/2 .* Prs(FB).^2 + (1-gamma)/2 .* Prs(TB).^2; 
PHI = (1+gamma)/2 .* (Prs_square(FB)-Prs_square(TB)) + (1-gamma)/2 .* (Prs_square(TB)-Prs_square(FB));
gasFlowSOCcons = [
    (gamma-1) .* gasFlowMax / 2 <= gasFlow_sum <= (gamma+1) .* gasFlowMax / 2;
    repmat((gamma-1) .* gasFlowMax / 2,[1,nGasType]) <= gasFlow <= repmat((gamma+1) .* gasFlowMax / 2, [1,nGasType]);
    Prsmin.^2 <= Prs_square <= Prsmax.^2;
    cone([PHI'+1; (2./newC .* gasFlow_sum)'; PHI'-1]);
    PHI <= (gasFlow_sum_ref.^2 + 2*gasFlow_sum_ref .* (gasFlow_sum - gasFlow_sum_ref)) ./ newC.^2 + sigma_PHI;
    sigma_PHI >= 0;
    ]:'gasFlowSOCcons';
% gas composition Taylor
PTGbus = mpc.ptg(:,1); 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbusMethane = CgsPTG * Qptg(:,1) ; QptgInbusHydrogen = CgsPTG * Qptg(:,2);
QptgInbusForAllGasComposition = [QptgInbusMethane, zeros(nGb,3),QptgInbusHydrogen,zeros(nGb,2)];
gasCompositionCons1 = [];

for r = 1:nGasType
    nodalGasInjectionForEachComp(:,r) = nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
            gamma, gasFlow(:,r),mpc,nGs,nGb,nGl);
end
gasCompositionCons = [
    -sigma_x <= repmat(W_ref,[1,nGasType]) .* gasComposition_ref + (repmat(sum(nodalGasInjectionForEachComp,2),[1,nGasType]) - repmat(W_ref,[1,nGasType])) .* gasComposition_ref ...
        + repmat(W_ref,[1,nGasType]) .* (gasComposition - gasComposition_ref) - nodalGasInjectionForEachComp <= sigma_x;
    - sigma_Gf <= (repmat((1+gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,1),:) + repmat((1-gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,2),:))/2 .* repmat(gasFlow_sum_ref,[1,nGasType]) ...
        + ( (repmat((1+gamma),[1,nGasType]).*gasComposition(mpc.Gline(:,1),:) + repmat((1-gamma),[1,nGasType]).*gasComposition(mpc.Gline(:,2),:))/2 - (repmat((1+gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,1),:) + repmat((1-gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,2),:))/2 ) .* repmat(gasFlow_sum_ref,[1,nGasType]) ...
        + (repmat((1+gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,1),:) + repmat((1-gamma),[1,nGasType]).*gasComposition_ref(mpc.Gline(:,2),:))/2 .* ( repmat(gasFlow_sum,[1,nGasType]) - repmat(gasFlow_sum_ref,[1,nGasType]) ) ...
        - gasFlow <= sigma_Gf;
    sigma_x >= 0;
    sum(gasComposition,2) == 1;
    gasComposition(12,:) == gasComposition(17,:);
    ]:'gasCompositionCons';
% summarize all the cons
constraints = [
    gasDemandCons;
    gasDemandBalanceCons;
    nodalGasFlowBalanceCons;
    PTGcons;
    GPPcons;
    electricityCons;
    electricityBalanceCons;
    WobbeIndexCons;
    FScons;
    SGcons;
    otherSecurityCons;
    gasFlowSOCcons;
    gasCompositionCons;
    gasSourceCons;
    Prs_square(20) == 25^2;
    ];
% 4. objective ============================================================
objfcn = obj_operatingCost(Pg,PGs,LCg,Qptg, mpc,CDF) ...
    +  1*alpha_PHI * sum(sum(sigma_PHI)) + 10000* alpha_x * sum(sum(sigma_x)) ...
    + alpha_Qd * 100 * (sum(sum(sigma_Qd)) + sum(sum(sigma_Qgpp)) + sum(sum(sigma_Gf)));
%%
% 要get dual，必须先设置好，而且把约束用cone函数写
options = sdpsettings('verbose',2,'solver','gurobi', 'debug',1,'usex0',0);
if solvedFlag == 0
    parameters = {gasComposition_ref,Qd_ref,Qgpp_ref,gasFlow_sum_ref,newC,W_ref,alpha_PHI,alpha_x,alpha_Qd};
    wantedVariables = {Prs_square, PGs, Qd, Qptg, Pptg, Pg, Pgpp, Qgpp, Va, gasComposition, gasFlow, gasFlow_sum, ...
        S_nodal, GCV_nodal, WI_nodal, FSnodal, PHI, sigma_PHI, sigma_x, sigma_Qd, sigma_Qgpp, sigma_Gf};
    GEopf_alternativeGas_optimizer = optimizer(constraints,objfcn,options,parameters,wantedVariables);
else
    GEopf_alternativeGas_optimizer = [];
end
end