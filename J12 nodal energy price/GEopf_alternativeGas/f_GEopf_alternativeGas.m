function [information, solution,solverTime] = f_GEopf_alternativeGas(mpc,options)
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
nGPP = size(mpc.gfuIndex,1);
for i = 1:nGPP
    iGppGd(i,1) = mpc.GEcon(find(mpc.GEcon(:,2) == mpc.gen(mpc.gfuIndex(i),1)),1);
end
nGs = size(mpc.Gsou,1);
nGd = size(find(mpc.Gbus(:,3)~=0),1);
nPTG = size(mpc.ptg,1);
iGd = find(mpc.Gbus(:,3)~=0);

% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
nGasType = 7; iCombustibleGas = 1:5; iNonCombustibleGas = 6:7;
[GCV, M, fs, a, R, T_stp, Prs_stp, Z_ref, T_gas, eta, CDF] = initializeParameters_J15();
%%
if  isempty(options)
    FSlimit = 0.2;
else
    FSlimit = options.FSlimit;
end


ob.Pg = [];

% step2: 初始状态认为是天然气的混合物
[GEinformation, GEsolution] = GEopf_milp(mpc);
if GEinformation.problem ~= 0
    error('optimization problem has issues');
end
% parameters for original natural gas
gasComposition_ng = mean(mpc.gasCompositionForGasSource);
M_ng = gasComposition_ng * M.all';
GCV_ng = gasComposition_ng * GCV.all';
S_ng = M_ng/M.air;
R_ng = gasComposition_ng * R.all';
WI_ng = GCV_ng / sqrt(S_ng);
FS_ng = gasComposition_ng * fs.All';
signGf = GEsolution.gamma;
%% form the optimization model
% parameters
gasComposition_ref = sdpvar(nGb,nGasType);
Qd_ref = sdpvar(nGd,nGasType);
Qgpp_ref = sdpvar(nGPP,nGasType);
gasFlow_sum_ref = sdpvar(nGl,1);
newC = sdpvar(nGl,1);
W_ref = sdpvar(nGb,1);
alpha_PHI = sdpvar(1,1);
alpha_x = sdpvar(1,1);
alpha_Qd = sdpvar(1,1);

% 5. optimizer ============================================================
solvedFlag = 0;
[~,~,GEopf_alternativeGas_optimizer] = GEopf_alternativeGas_construct(gasComposition_ref, Qd_ref, Qgpp_ref, gasFlow_sum_ref, newC, W_ref, alpha_PHI, alpha_x, alpha_Qd, ...
    mpc,nb,nGb,nGl,nGs,nGd,nGasType,nPTG,nGPP,iGd,iGppGd,signGf,REF,GCV,M,S_ng,FSlimit,FS_ng,GCV_ng,eta, ...
    fs,CDF,solvedFlag);


%% initial value
% 需要提前计算并迭代的变量：
PGs_ref = GEsolution.PGs;
S_pipeline_ref = S_ng; %
R_pipeline_ref = R_ng;
gasComposition_ref = repmat(gasComposition_ng,[nGb,1]);
Qd_ref = repmat(mpc.Gbus(iGd,3),[1,nGasType]) .* gasComposition_ref(iGd,:);
Qgpp_ref = repmat(GEsolution.Qgpp,[1,nGasType]) .* gasComposition_ref(iGppGd,:);
gasFlow_sum_ref = GEsolution.Gf;
newC = mpc.Gline(:,3) .* sqrt(R_ng) ./ sqrt(R_pipeline_ref);
PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbus_ref = CgsPTG * GEsolution.Qptg; 
W_ref = nodalGasInjection(PGs_ref,QptgInbus_ref,signGf, gasFlow_sum_ref,mpc,nGs,nGb,nGl);
gasLoad_ng = mpc.Gbus(iGd,3);

%% sequential programming
alpha_PHI = 1e-1;alpha_x = 1e-1; alpha_Qd = 1e-1;
lambda = 1; % coefficient for gas flow
multiplier = 1.5; alpha_PHI_max = 1e4; alpha_x_max = 1e4; alpha_Qd_max = 1e4;
iterationMax = 60;
gap = 1e-3; % optimization convergence criterion
solverTime = 0;
for v = 1:iterationMax
tic
[xvalue, errorcode{v}] = GEopf_alternativeGas_optimizer(gasComposition_ref,Qd_ref,Qgpp_ref,gasFlow_sum_ref,newC,W_ref,alpha_PHI,alpha_x,alpha_Qd);
solverTime = solverTime + toc;
[Prs_square, PGs, Qd, Qptg, Pptg, Pg, Pgpp, Qgpp, Va, gasComposition, gasFlow, gasFlow_sum, ...
    S_nodal, GCV_nodal, WI_nodal, FSnodal, PHI, sigma_PHI, sigma_x, sigma_Qd, sigma_Qgpp, sigma_Gf] = deal(xvalue{:});
Prs = sqrt(Prs_square);
%% results
LCg = zeros(nGd,1);gamma = signGf;
[objfcn,totalCost,electricityGenerationCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg, mpc, CDF);
%
[sol{v}.totalCost,sol{v}.genAndLCeCost,sol{v}.gasPurchasingCost,sol{v}.gasCurtailmentCost,sol{v}.PTGsubsidy,sol{v}.objfcn, ...
    sol{v}.Prs_square,sol{v}.Prs,sol{v}.PGs,sol{v}.Qd,sol{v}.Qptg,sol{v}.Pptg,sol{v}.Pg,sol{v}.Pgpp,sol{v}.Qgpp,...
    sol{v}.Va,sol{v}.LCg,sol{v}.gamma,sol{v}.gasComposition,sol{v}.gasFlow,sol{v}.gasFlow_sum,sol{v}.S_nodal,...
    sol{v}.GCV_nodal,sol{v}.WI,sol{v}.PHI,sol{v}.sigma_PHI,sol{v}.sigma_x,sol{v}.sigma_Qd,sol{v}.sigma_Qgpp,sol{v}.FSnodal] = ...
    deal(totalCost,electricityGenerationCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy, objfcn, ...
    Prs_square,Prs,PGs,Qd,Qptg,Pptg,Pg,Pgpp,Qgpp,Va,LCg,gamma,gasComposition,gasFlow,gasFlow_sum,S_nodal,...
    GCV_nodal,WI_nodal,PHI,sigma_PHI,sigma_x,sigma_Qd,sigma_Qgpp,FSnodal);

%% covergence criterion
if v > 1
    criterion.sigma_PHI(v-1) = sum(sum(sigma_PHI))/1e2;
    criterion.sigma_x(v-1) = sum(sum(sigma_x));
    criterion.sigma_Qd(v-1) = sum(sum(sigma_Qd));
    criterion.sigma_Qgpp(v-1) = sum(sum(sigma_Qgpp));
    criterion.delta_sigma_PHI(v-1) = abs( (sum(sum(sigma_PHI)) - sum(sum(sol{v-1}.sigma_PHI)))) / abs( (sum(sum(sigma_PHI)) + sum(sum(sol{v-1}.sigma_PHI))));
    criterion.delta_sigma_x(v-1) = abs( (sum(sum(sigma_x)) - sum(sum(sol{v-1}.sigma_x)))) / abs( (sum(sum(sigma_x)) + sum(sum(sol{v-1}.sigma_x))));
    criterion.delta_sigma_Qd(v-1) = abs( (sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)))) / abs( (sum(sum(sigma_Qd)) + sum(sum(sol{v-1}.sigma_Qd))));
    criterion.delta_sigma_Qgpp(v-1) = abs( (sum(sum(sigma_Qgpp)) - sum(sum(sol{v-1}.sigma_Qgpp)))) / abs( (sum(sum(sigma_Qgpp)) + sum(sum(sol{v-1}.sigma_Qgpp))));
    %cost
    criterion.totalCost(v-1) = abs( (totalCost - sol{v-1}.totalCost)) / sol{v-1}.totalCost;
    criterion.genAndLCeCost(v-1) = abs( (electricityGenerationCost - sol{v-1}.genAndLCeCost)/sol{v-1}.genAndLCeCost );
    criterion.gasPurchasingCost(v-1) = abs( (gasPurchasingCost - sol{v-1}.gasPurchasingCost)/sol{v-1}.gasPurchasingCost );
    criterion.PTGsubsidy(v-1) = abs( (PTGsubsidy - sol{v-1}.PTGsubsidy)/sol{v-1}.PTGsubsidy );
    criterion.gasCompositionCH4(v-1) = abs(sum(gasComposition(:,1)) - sum(sol{v-1}.gasComposition(:,1)))*2 ./ (sum(gasComposition(:,1)) + sum(sol{v-1}.gasComposition(:,1)));
    criterion.gasCompositionH2(v-1) = abs(sum(gasComposition(:,5)) - sum(sol{v-1}.gasComposition(:,5)))*2 / (sum(gasComposition(:,5)) + sum(sol{v-1}.gasComposition(:,5)));
    criterion.S_nodal(v-1) = max(abs( (S_nodal-sol{v-1}.S_nodal)./sol{v-1}.S_nodal ));
    criterion.Qd(v-1) = abs( ( sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)) ) / sum(sum(sol{v-1}.sigma_Qd)) );

    if  ( ((criterion.sigma_PHI(v-1)<=1e-3) && (criterion.sigma_x(v-1)<=1e-2) && (criterion.sigma_Qd(v-1)<=1e-2) && (criterion.sigma_Qgpp(v-1)<=1e-2) ) && ...
                ( (criterion.sigma_x(v-1)<=gap) || (criterion.delta_sigma_x(v-1)<=gap) ) ...
                && ( (criterion.sigma_Qd(v-1)<=gap) ||  (criterion.delta_sigma_Qd(v-1)<gap) ) ...
                && ( (criterion.sigma_Qgpp(v-1)<=gap) || (criterion.delta_sigma_Qgpp(v-1)<gap) ) ) || v == iterationMax
            % get dual
            yalmip('clear')
            solvedFlag = 1;
            [objfcn,constraints,~] = GEopf_alternativeGas_construct(gasComposition_ref, Qd_ref, Qgpp_ref, gasFlow_sum_ref, newC, W_ref, alpha_PHI, alpha_x, alpha_Qd, ...
                mpc,nb,nGb,nGl,nGs,nGd,nGasType,nPTG,nGPP,iGd,iGppGd,signGf,REF,GCV,M,S_ng,FSlimit,FS_ng,GCV_ng,eta, ...
                fs,CDF,solvedFlag);
            options_new = sdpsettings('verbose',2,'solver','mosek', 'debug',1);
            optimize(constraints,objfcn,options_new);
            nodalGasPriceForEachComposition = - dual(constraints('nodalGasFlowBalanceCons'));
            nodalGasPrice = sum(nodalGasPriceForEachComposition .* gasComposition,2); % $/Mm3
            nodalElectricityPrice = dual(constraints('electricityBalanceCons'))/100;
            sol{v}.nodalGasPriceForEachComposition = nodalGasPriceForEachComposition;
            sol{v}.nodalGasPrice = nodalGasPrice;
            sol{v}.nodalElectricityPrice = nodalElectricityPrice;
            break
        else
        alpha_PHI = min([multiplier * alpha_PHI,alpha_PHI_max]);
        alpha_x = min([multiplier * alpha_x,alpha_x_max]);
        alpha_Qd = min([multiplier * alpha_Qd,alpha_Qd_max]);
    end
end
%% update S and Z
R_nodal = gasComposition * R.all';
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
R_pipeline_ref = ( (1+gamma).*R_nodal(FB) + (1-gamma).*R_nodal(TB) ) / 2;
QptgInbus = CgsPTG * sum(Qptg,2);
W_ref = nodalGasInjection(PGs,QptgInbus,gamma,gasFlow_sum,mpc,nGs,nGb,nGl);
gasComposition_ref = gasComposition;
gasFlow_sum_ref = gasFlow_sum;
newC = mpc.Gline(:,3) .* sqrt(R_ng) ./ sqrt(R_pipeline_ref);

end
information.errorcode = errorcode{v}; solution = sol{v};
if v==iterationMax % if do not converge
    information.ruleoutFlag = 1;
else
    information.ruleoutFlag = 0;
end







