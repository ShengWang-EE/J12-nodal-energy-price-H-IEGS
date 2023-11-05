function [information, solution] = GEopf_milp(mpc)
% 20220810: 纯电力天然气的稳态opf，考虑燃气机组和ptg，考虑切负荷
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
% create (read-only) copies of individual fields for convenience\
baseMVA=100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
%% initialization
    gasLoadCurtailmentProportion = 0;
    electricityLoadCurtailmentProportion = 0;

nGb = size(mpc.Gbus,1);
nGs = size(mpc.Gsou,1);
nGl = size(mpc.Gline,1);
iGd = find(mpc.Gbus(:,3)~=0);
nGd = size(iGd,1);
nb = size(mpc.bus,1);
ng = size(mpc.gen,1);
id = find(mpc.bus(:,3)~=0);
nd = size(id,1);
nPTG = size(mpc.ptg,1);

refbus = find(mpc.bus(:, BUS_TYPE) == REF);

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF] = initializeParameters_J15();
%% state variable
Prs_square = sdpvar(nGb,1); % nodal gas pressure
PGs = sdpvar(nGs,1); % gas production of gas source
Gf = sdpvar(nGl,1); % gas flow in the pipeline
LCg = sdpvar(nGd,1); % nodal gas load curtailment
Qptg = sdpvar(nPTG,1); % gas (hydrogen) production of PTG

Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)
LCe = sdpvar(nd,1); % electricity load curtailment

gamma = binvar(nGl,1);% direction of gas flow, 0-1
gamma = (gamma-0.5)*2; % 1,-1
PHI = sdpvar(nGl,1); %auxiliary variable for gas flow
%% upper and lower bounds
Prs_square_min = mpc.Gbus(:,5).^2; Prs_square_max = mpc.Gbus(:,6).^2;
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4);
Gfmin = -mpc.Gline(:,5); Gfmax = mpc.Gline(:,5); 
LCgmin = 0; LCgmax = gasLoadCurtailmentProportion * mpc.Gbus(iGd,3);
Qptgmin = 0; Qptgmax = mpc.ptg(:,6);
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
LCemin = 0; LCemax = electricityLoadCurtailmentProportion * mpc.bus(id,PD)/baseMVA;

% Qptg = 0.5*Qptgmax;
%% contraints
Pptg = Qptg/24/3600 * GCV.ng_ref * eta.electrolysis / baseMVA;
boxCons = [
    Prs_square_min <= Prs_square <= Prs_square_max;
    PGsmin <= PGs <= PGsmax;
    LCgmin <= LCg <= LCgmax;
    Qptgmin <= Qptg <= Qptgmax;
%     Vamin <= Va <= Vamax;
    Pgmin <= Pg <= Pgmax;
    LCemin <= LCe <= LCemax;
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,LCe,mpc,id) == 0;];
electricityBranchFlowConsDC = [consfcn_electricityBranchFlow(Va, mpc, il) <= 0;];
Pgpp = Pg(mpc.gfuIndex)*baseMVA;
Qgpp = Pgpp / GCV.ng_ref * 3600 * 24 * eta.GFU;
gasNodalBalanceCons = [consfcn_gasNodalBalance(Gf,Pg,PGs,Qgpp,LCg,Qptg,mpc,iGd,nGb,nGl,nGs,nGd) == 0;];

gasPipelineFlowCons = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];
% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
gasFlowCons = [
    Prs_square(FB) - Prs_square(TB) == Gf.^2 ./ mpc.Gline(:,3).^2;
    ];

%% solotion
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasPipelineFlowCons;
    gasFlowCons;
    Prs_square(20) == 25^2;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,LCe,PGs,Qptg,LCg,mpc,CDF.electricity,CDF.gas,id,iGd);
yalmipOptions = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1);

information = optimize(cons,objfcn,yalmipOptions);
%% results
Prs = sqrt(value(Prs_square)); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
LCg = value(LCg); % nodal gas load curtailment
Qptg = value(Qptg); % gas (hydrogen) production of PTG
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
LCe = value(LCe); % electricity load curtailment
gamma = value(gamma);
objfcn = value(objfcn);
Qgpp = value(Qgpp);
[operatingCost,electricityGenerationCost,electricityLoadCurtailmentCost,...
    gasPurchasingCost,gasCurtailmentCost] = objfcn_IEGSoperatingCost(Pg,LCe,PGs,Qptg,LCg,mpc,CDF.electricity,CDF.gas,id,iGd);

[solution.Prs, solution.PGs, solution.Gf, solution.LCg, solution.Qptg, ...
    solution.Va, solution.Pg, solution.LCe, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost, solution.electricityLoadCurtailmentCost, ...
    solution.gasPurchasingCost, solution.gasCurtailmentCost] = deal(...
    Prs, PGs, Gf, LCg, Qptg, Va, Pg, LCe, gamma, ...
    operatingCost, electricityGenerationCost, electricityLoadCurtailmentCost, ...
    gasPurchasingCost, gasCurtailmentCost);
solution.Qgpp = Qgpp;
end
