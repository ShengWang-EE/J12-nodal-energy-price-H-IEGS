function [information, solution] = GEopf_alternativeGas_optimizer(mpc)
%
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
%
nb   = size(mpc.bus, 1);    %% number of buses
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl = size(mpc.Gline,1);
nGen = sum(mpc.gen,1);
nGpp = size(mpc.gfuIndex,1);
for i = 1:nGpp
    iGppGd(i,1) = mpc.GEcon(find(mpc.GEcon(:,2) == mpc.gen(mpc.gfuIndex(i),1)),1);
end
nGs = size(mpc.Gsou,1);
nGasLoad = size(find(mpc.Gbus(:,3)~=0),1);
nPTG = size(mpc.ptg,1);
iGd = find(mpc.Gbus(:,3)~=0);

% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
nGasType = 7; iCombustibleGas = 1:5; iNonCombustibleGas = 6:7;
[GCV, M, fs, a, R, T_stp, Prs_stp, Z_ref, T_gas, eta, CDF] = initializeParameters_J15();
%%
xi = 0.1; % threshold for the security index
alpha_PHI = 1e-4;alpha_x = 1e-4; alpha_Qd = 1e-4;
lambda = 1; % coefficient for gas flow
multiplier = 1.5; alpha_PHI_max = 1e2; alpha_x_max = 1e4; alpha_Qd_max = 1e4;
iterationMax = 60;
gap = 1e-3; % optimization convergence criterion

ob.Pg = [];

% step2: 初始状态认为是天然气的混合物
[GEinformation, GEsolution] = GEopf_milp(mpc);
if GEinformation.problem ~= 0
    error('optimization problem has issues');
end
%% initial value
PGs_ref = GEsolution.PGs;
% parameters for original natural gas
gasComposition_ng = mean(mpc.gasCompositionForGasSource);
M_ng = gasComposition_ng * M.all';
GCV_ng = gasComposition_ng * GCV.all';
S_ng = M_ng/M.air;
R_ng = gasComposition_ng * R.all';
WI_ng = GCV_ng / sqrt(S_ng);
FS_ng = gasComposition_ng * fs.All';
% 需要提前计算并迭代的变量：
S_pipeline_ref = S_ng; %
R_pipeline_ref = R_ng;
gasComposition_ref = repmat(gasComposition_ng,[nGb,1]);
Qd_ref = repmat(mpc.Gbus(iGd,3),[1,nGasType]) .* gasComposition_ref(iGd,:);
Qgpp_ref = repmat(GEsolution.Qgpp,[1,nGasType]) .* gasComposition_ref(iGppGd,:);
gasFlow_sum_ref = GEsolution.Gf;

PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbus_ref = CgsPTG * GEsolution.Qptg; 

signGf = GEsolution.gamma;
W_ref = nodalGasInjection(PGs_ref,QptgInbus_ref,signGf, gasFlow_sum_ref,mpc,nGs,nGb,nGl);
gasLoad_ng = mpc.Gbus(iGd,3);

solverTime = zeros(iterationMax,1);

%% optimizer

%%
for v = 1:iterationMax

%% state variables
Prs_square = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGasLoad,nGasType);% Mm3/day
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = sdpvar(nPTG,1); % electricity consumption, 1/100 MW
Pg = sdpvar(size(mpc.gen,1),1); % include TPP, GPP and renewable generators, 1/100 MW
Qgpp = sdpvar(nGpp, nGasType); % Mm3/day
Va = sdpvar(nb,1);
gasComposition = sdpvar(nGb,nGasType); 
gasFlow = sdpvar(nGl,nGasType);% Mm3/day
PHI = sdpvar(nGl,1); %auxiliary variable for gas flow
sigma_PHI = sdpvar(nGl,1); % error limit for gas flow
sigma_x = sdpvar(nGb, nGasType); % error for gas composition
sigma_Qd = sdpvar(nGasLoad, nGasType);
sigma_Qgpp = sdpvar(nGpp,nGasType);
% ------------test------------------
LCg = zeros(nGasLoad,1);
% varepsilon_x = zeros(nGb, nGasType); % error for gas composition
gamma = signGf;
%% bounds
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6); % bar
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4); % Mm3/day
Qptgmin = mpc.ptg(:,5); Qptgmax = mpc.ptg(:,6);
QptgMax_hydrogen = Qptgmax;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
Pgmax(34:end,:) = 1e-4; % no electricity load curtailment
LCgmin = zeros(nGasLoad,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0;  
refs = find(mpc.bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = 1;   %% voltage angle reference constraints
Val(refs) = 1;
gasFlowMax = mpc.Gline(:,5);
%% pre-calculation
% calculate the gas flow of each pipeline
newC = mpc.Gline(:,3) .* sqrt(R_ng) ./ sqrt(R_pipeline_ref);
% newC = mpc.Gline(:,3);
%% constraints
% gas source cons
gasSourceCons = [  PGsmin <= PGs <= PGsmax;];
% gas demand cons
energyDemand = mpc.Gbus(iGd,3) * GCV_ng; % energy need of these gas bus
gasDemandCons = [
     -sigma_Qd <= gasComposition_ref(iGd,:) .* repmat(sum(Qd_ref,2),[1,nGasType]) + (gasComposition(iGd,:)-gasComposition_ref(iGd,:)) .* repmat(sum(Qd_ref,2),[1,nGasType]) ...
        + gasComposition_ref(iGd,:) .* (repmat(sum(Qd,2),[1,nGasType])-repmat(sum(Qd_ref,2),[1,nGasType])) - Qd <= sigma_Qd;
    Qd * GCV.all'/1e9 == energyDemand/1e9;    
    Qd >= 0;
    sigma_Qd >= 0;
    ]:'gasDemandCons';
% nodal gas flow balance cons
nodalGasFlowBalanceCons = [
    consfcn_nodalGasFlowBalance(PGs,Qd,Qgpp,Qptg, gasFlow,mpc,nGasType,nGpp,nGasLoad,iGd) == 0;
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
    Pgpp == Qgpp/24/3600 * GCV.all';
    Qgpp >= 0;
    ]:'GPPcons';
% electricity flow
electricityCons = [
    consfcn_electricPowerBalance(Va,Pg,Pptg,mpc) == 0;
    consfcn_electricBranchFlow(Va, mpc, il) <= 0;
    Pgmin <= Pg <= Pgmax;
    Va(refs) == 0; % 除了slackbus外，其他相角都没约束。但是一些solver在处理inf的上下限的时候有问题
    ]:'electricityCons';
% wobbe index
GCV_nodal = gasComposition * GCV.all';
S_nodal = gasComposition * M.all' / M.air;
sqrtS = 0.5 * (S_nodal/sqrt(S_ng) + sqrt(S_ng));
WI_nodal = GCV_nodal ./ sqrtS; % 如果不能自动转化，那就手动化一下
% WImin = 47.2; WImax = 51.41; %old
WImin = 46.5; WImax = 52.85;
WobbeIndexCons = [
    WImin * sqrtS <= GCV_nodal/1e6 <= WImax * sqrtS
    ]:'WobbeIndexCons';
% Weaver flame speed factor
FSnodal = gasComposition * fs.All';
FScons = [
    (1-xi) * FS_ng <= FSnodal <= (1+xi) * FS_ng;
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
PHI = (1+gamma)/2 .* (Prs_square(FB)-Prs_square(TB)) + (1-gamma)/2 .* (Prs_square(TB)-Prs_square(FB));
gasFlowSOCcons = [
    (gamma-1) .* gasFlowMax / 2 <= gasFlow_sum <= (gamma+1) .* gasFlowMax / 2;
    repmat((gamma-1) .* gasFlowMax / 2,[1,nGasType]) <= gasFlow <= repmat((gamma+1) .* gasFlowMax / 2, [1,nGasType]);
    Prsmin.^2 <= Prs_square <= Prsmax.^2;
    PHI >= gasFlow_sum.^2 ./ newC.^2 ;
    PHI <= (gasFlow_sum_ref.^2 + 2*gasFlow_sum_ref .* (gasFlow_sum - gasFlow_sum_ref)) ./ newC.^2 + sigma_PHI;
    sigma_PHI >= 0;
    ]:'gasFlowSOCcons';
% gas composition Taylor
PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbusMethane = CgsPTG * Qptg(:,1) ; QptgInbusHydrogen = CgsPTG * Qptg(:,2);
QptgInbusForAllGasComposition = [QptgInbusMethane, zeros(nGb,3),QptgInbusHydrogen,zeros(nGb,2)];
gasCompositionCons1 = [];

for r = 1:nGasType
    nodalGasInjectionForEachComp(:,r) = nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
            gamma, gasFlow(:,r),mpc,nGs,nGb,nGl);
end

% for r = 1:nGasType
%     gasCompositionCons1 = [
%         gasCompositionCons1;
% %         gasComposition(:,r) == ...
% %             nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
% %             gamma, gasFlow(:,r),mpc,nGs,nGb,nGl) ./ W_ref + varepsilon_x(:,r);
% %         gasComposition(:,r) == nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
% %             gamma, gasFlow(:,r),mpc,nGs,nGb,nGl) ./ nodalGasInjectionSum; % original cons
%         ];
% end
gasCompositionCons = [
    -sigma_x <= repmat(W_ref,[1,nGasType]) .* gasComposition_ref + (repmat(sum(nodalGasInjectionForEachComp,2),[1,nGasType]) - repmat(W_ref,[1,nGasType])) .* gasComposition_ref ...
        + repmat(W_ref,[1,nGasType]) .* (gasComposition - gasComposition_ref) - nodalGasInjectionForEachComp <= sigma_x;
    sigma_x >= 0;
    sum(gasComposition,2) == 1;
    gasComposition(12,:) == gasComposition(17,:);
    ]:'gasCompositionCons';
% summarize all the cons
constraints = [
    gasDemandCons;
    nodalGasFlowBalanceCons;
    PTGcons;
    GPPcons;
    electricityCons;
    WobbeIndexCons;
    FScons;
    SGcons;
    otherSecurityCons;
    gasFlowSOCcons;
    gasCompositionCons;
    gasSourceCons;
    Prs_square(20) == 25^2;
    ];
%% solve the problem
objfcn = obj_operatingCost(Pg,PGs,LCg,Qptg,PHI, sigma_PHI, sigma_x, alpha_PHI,alpha_x, mpc,lambda,CDF) ...
    + alpha_Qd * 100 * (sum(sum(sigma_Qd))+sum(sum(sigma_Qgpp)));
options = sdpsettings('verbose',2,'solver','gurobi', 'debug',1,'usex0',0);
% options.ipopt.tol = 1e-4;
output{v} = optimize(constraints, objfcn, options);

solverTime(v) = output{v}.solvertime;
%% results
Prs_square = value(Prs_square);
Prs = sqrt(value(Prs_square));
PGs = value(PGs); % Mm3/day
Qd = value(Qd); % Mm3/day
Qptg = value(Qptg);
Pptg = value(Pptg) * baseMVA; % MW
Pg = value(Pg) * baseMVA;  % MW
Pgpp = value(Pgpp) * baseMVA;  % MW
Qgpp = value(Qgpp);  % Mm3/day
Va = value(Va);
LCg = value(LCg);
gamma = value(gamma);
gasComposition = value(gasComposition); % hydrogen, gas
gasFlow = value(gasFlow);
gasFlow_sum = sum(value(gasFlow_sum),2);
S_nodal = value(S_nodal);
GCV_nodal = value(GCV_nodal)/1e6; % MJ/m3
WInodal = value(WI_nodal)/1e6; % MJ/m3
FSnodal = value(FSnodal);
PHI = value(PHI);
sigma_PHI = value(sigma_PHI);
sigma_x = value(sigma_x);
sigma_Qd = value(sigma_Qd);
sigma_Qgpp = value(sigma_Qgpp);

[objfcn,totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy,penalty_PHI,penalty_sigma_PHI, penalty_sigma_x] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg,PHI, sigma_PHI, sigma_x, alpha_PHI,alpha_x, mpc,lambda, CDF);
%
[sol{v}.totalCost,sol{v}.genAndLCeCost,sol{v}.gasPurchasingCost,sol{v}.gasCurtailmentCost,sol{v}.PTGsubsidy,sol{v}.penalty_PHI,sol{v}.penalty_sinma_PHI, sol{v}.penalty_sigma_x, sol{v}.objfcn, ...
    sol{v}.Prs_square,sol{v}.Prs,sol{v}.PGs,sol{v}.Qd,sol{v}.Qptg,sol{v}.Pptg,sol{v}.Pg,sol{v}.Pgpp,sol{v}.Qgpp,...
    sol{v}.Va,sol{v}.LCg,sol{v}.gamma,sol{v}.gasComposition,sol{v}.gasFlow,sol{v}.gasFlow_sum,sol{v}.S_nodal,...
    sol{v}.GCV_nodal,sol{v}.WI,sol{v}.PHI,sol{v}.sigma_PHI,sol{v}.sigma_x,sol{v}.sigma_Qd,sol{v}.sigma_Qgpp] = ...
    deal(totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy,penalty_PHI,penalty_sigma_PHI, penalty_sigma_x, objfcn, ...
    Prs_square,Prs,PGs,Qd,Qptg,Pptg,Pg,Pgpp,Qgpp,Va,LCg,gamma,gasComposition,gasFlow,gasFlow_sum,S_nodal,...
    GCV_nodal,WInodal,PHI,sigma_PHI,sigma_x,sigma_Qd,sigma_Qgpp);

%% covergence criterion
if v > 1
    criterion.sigma_PHI(v-1) = sum(sum(sigma_PHI))/1e3;
    criterion.sigma_x(v-1) = sum(sum(sigma_x));
    criterion.sigma_Qd(v-1) = sum(sum(sigma_Qd));
    criterion.sigma_Qgpp(v-1) = sum(sum(sigma_Qgpp));
    criterion.delta_sigma_PHI(v-1) = abs( (sum(sum(sigma_PHI)) - sum(sum(sol{v-1}.sigma_PHI)))) / abs( (sum(sum(sigma_PHI)) + sum(sum(sol{v-1}.sigma_PHI))));
    criterion.delta_sigma_x(v-1) = abs( (sum(sum(sigma_x)) - sum(sum(sol{v-1}.sigma_x)))) / abs( (sum(sum(sigma_x)) + sum(sum(sol{v-1}.sigma_x))));
    criterion.delta_sigma_Qd(v-1) = abs( (sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)))) / abs( (sum(sum(sigma_Qd)) + sum(sum(sol{v-1}.sigma_Qd))));
    criterion.delta_sigma_Qgpp(v-1) = abs( (sum(sum(sigma_Qgpp)) - sum(sum(sol{v-1}.sigma_Qgpp)))) / abs( (sum(sum(sigma_Qgpp)) + sum(sum(sol{v-1}.sigma_Qgpp))));
    %cost
    criterion.totalCost(v-1) = abs( (totalCost - sol{v-1}.totalCost)) / sol{v-1}.totalCost;
    criterion.genAndLCeCost(v-1) = abs( (genAndLCeCost - sol{v-1}.genAndLCeCost)/sol{v-1}.genAndLCeCost );
    criterion.gasPurchasingCost(v-1) = abs( (gasPurchasingCost - sol{v-1}.gasPurchasingCost)/sol{v-1}.gasPurchasingCost );
    criterion.PTGsubsidy(v-1) = abs( (PTGsubsidy - sol{v-1}.PTGsubsidy)/sol{v-1}.PTGsubsidy );
    criterion.gasCompositionCH4(v-1) = abs(sum(gasComposition(:,1)) - sum(sol{v-1}.gasComposition(:,1)))*2 ./ (sum(gasComposition(:,1)) + sum(sol{v-1}.gasComposition(:,1)));
    criterion.gasCompositionH2(v-1) = abs(sum(gasComposition(:,5)) - sum(sol{v-1}.gasComposition(:,5)))*2 / (sum(gasComposition(:,5)) + sum(sol{v-1}.gasComposition(:,5)));
    criterion.S_nodal(v-1) = max(abs( (S_nodal-sol{v-1}.S_nodal)./sol{v-1}.S_nodal ));
    criterion.Qd(v-1) = abs( ( sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)) ) / sum(sum(sol{v-1}.sigma_Qd)) );

    if ( (criterion.sigma_PHI(v-1)<=gap) || (criterion.delta_sigma_PHI(v-1)<=gap) ) ...
            && ( (criterion.sigma_x(v-1)<=gap) || (criterion.delta_sigma_x(v-1)<=gap) ) ...
            && ( (criterion.sigma_Qd(v-1)<=gap) ||  (criterion.delta_sigma_Qd(v-1)<gap) ) ...
            && ( (criterion.sigma_Qgpp(v-1)<=gap) || (criterion.delta_sigma_Qgpp(v-1)<gap) )
        break
    else
        alpha_PHI = min([multiplier * alpha_PHI,alpha_PHI_max]);
        alpha_x = min([multiplier * alpha_x,alpha_x_max]);
        alpha_Qd = min([multiplier * alpha_Qd,alpha_Qd_max]);
    end
end
% ob
ob.Pg = [ob.Pg, Pg];
%% update S and Z
R_nodal = gasComposition * R.all';
R_pipeline_ref = ( (1+gamma).*R_nodal(FB) + (1-gamma).*R_nodal(TB) ) / 2;
QptgInbus = CgsPTG * sum(Qptg,2);
W_ref = nodalGasInjection(PGs,QptgInbus,gamma,gasFlow_sum,mpc,nGs,nGb,nGl);
gasComposition_ref = gasComposition;
if output{v}.problem ~= 0
%     error('optimization failed');
%     break
end
end
information = output{v}; solution = sol{v};
end







