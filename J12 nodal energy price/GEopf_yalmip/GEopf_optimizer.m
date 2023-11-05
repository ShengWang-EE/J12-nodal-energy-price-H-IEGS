function [results, information, solution] = GEopf_optimizer(mpc,GEoptions)
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
if isempty(GEoptions)
    GEoptions.gasLoadCurtailmentProportion = 1;
    GEoptions.electricityLoadCurtailmentProportion = 1;
end

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
penalty = 1e6;% penalty factor for SOC relaxations

eta_electrolysis = 0.7;    % from the energy perspective, the effciency is about 80%
eta_methanation = 0.8;
etaGFU = 0.4211;           % from the energy perspective, 从1/200换算而来
GCVng = 41.04 * 1e6;     % J/m3
CDFe = 1e4; % MW/hour, 大概数值，从jia文章中拿的
CDFg = CDFe * GCVng / 3600 / 24;
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
LCgmin = 0; LCgmax = GEoptions.gasLoadCurtailmentProportion * mpc.Gbus(iGd,3);
Qptgmin = 0; Qptgmax = mpc.ptg(:,6);
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
LCemin = 0; LCemax = GEoptions.electricityLoadCurtailmentProportion * mpc.bus(id,PD)/baseMVA;
%% contraints
Pptg = Qptg/24/3600 * GCVng * eta_electrolysis / baseMVA;
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
gasNodalBalanceCons = [consfcn_gasNodalBalance(Gf,Pg,PGs,LCg,Qptg,mpc,iGd,nGb,nGl,nGs,nGd,GCVng) == 0;];

gasPipelineFlowCons = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];
% SOC for gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
gasFlowCons = [
    PHI >= Gf.^2 ./ mpc.Gline(:,3).^2 ;
    PHI >= Prs_square(TB) - Prs_square(FB) + (gamma + 1) .* (Prs_square_min(FB).^2 - Prs_square_max(TB).^2);
    PHI >= Prs_square(FB) - Prs_square(TB) + (gamma - 1) .* (Prs_square_max(FB).^2 - Prs_square_min(TB).^2);
    PHI <= Prs_square(TB) - Prs_square(FB) + (gamma + 1) .* (Prs_square_max(FB).^2 - Prs_square_min(TB).^2);
    PHI <= Prs_square(FB) - Prs_square(TB) + (gamma - 1) .* (Prs_square_min(FB).^2 - Prs_square_max(TB).^2);
    ];

%% solotion
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasPipelineFlowCons;
    gasFlowCons;
    ];
penaltyTerm = penalty * sum(PHI);
objfcn = objfcn_IEGSoperatingCost(Pg,LCe,PGs,LCg,mpc,CDFe,CDFg) + penaltyTerm;
yalmipOptions = sdpsettings('verbose',2,'solver','mosek','usex0',0,'debug',1);

information = optimize(cons,objfcn,yalmipOptions);
%% results
results = mpc;

Prs = sqrt(value(Prs_square)); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
LCg = value(LCg); % nodal gas load curtailment
Qptg = value(Qptg); % gas (hydrogen) production of PTG
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
LCe = value(LCe); % electricity load curtailment
gamma = value(gamma);
[objfcn,electricityGenerationCost,electricityLoadCurtailmentCost,...
    gasPurchasingCost,gasCurtailmentCost] = objfcn_IEGSoperatingCost(Pg,LCe,PGs,LCg,mpc,CDFe,CDFg);

[solution.Prs, solution.PGs, solution.Gf, solution.LCg, solution.Qptg, ...
    solution.Va, solution.Pg, solution.LCe, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost, solution.electricityLoadCurtailmentCost, ...
    solution.gasPurchasingCost, solution.gasCurtailmentCost,solution.penaltyTerm] = deal(...
    Prs, PGs, Gf, LCg, Qptg, Va, Pg, LCe, gamma, ...
    objfcn, electricityGenerationCost, electricityLoadCurtailmentCost, ...
    gasPurchasingCost, gasCurtailmentCost,penaltyTerm);

results.Gbus(:,7) = Prs;
results.Gsou(:,5) = PGs;
results.Gline(:,6) = Gf;
results.Gbus(iGd,10) = LCg;
results.ptg(:,7) = Qptg;
results.bus(:,VA) = Va;
results.gen(:,PG) = Pg * baseMVA;
results.bus(id,14) = LCe;
[results.objfcn, results.electricityGenerationCost, results.electricityLoadCurtailmentCost, ...
    results.gasPurchasingCost, results.gasCurtailmentCost,results.penaltyTerm] = deal(...
    objfcn, electricityGenerationCost, electricityLoadCurtailmentCost, ...
    gasPurchasingCost, gasCurtailmentCost,penaltyTerm);

end
