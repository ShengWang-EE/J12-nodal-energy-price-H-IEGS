function [mpc, gtd] = case24GE_J15()
% v6:该案例对应的是电氢联合优化潮流，有ptg和windfarm
% 并且在mpce基础上修改ptg位置，来人为制造阻塞
% v7: 换了ptg的接入点
% 20221230：删去了dispatchable load，还是写在优化模型里作为单独的变量比较好，否则容易造成误解

%CASE24_IEEE_RTS  Power flow data for the IEEE RELIABILITY TEST SYSTEM.
%   Please see CASEFORMAT for details on the case file format.
%
%   This system data is from the IEEE RELIABILITY TEST SYSTEM, see
%
%   IEEE Reliability Test System Task Force of the Applications of
%   Probability Methods Subcommittee, "IEEE reliability test system,"
%   IEEE Transactions on Power Apparatus and Systems, Vol. 98, No. 6,
%   Nov./Dec. 1979, pp. 2047-2054.
%
%   IEEE Reliability Test System Task Force of Applications of
%   Probability Methods Subcommittee, "IEEE reliability test system-96,"
%   IEEE Transactions on Power Systems, Vol. 14, No. 3, Aug. 1999,
%   pp. 1010-1020.
%
%   Cost data is from Web site run by Georgia Tech Power Systems Control
%   and Automation Laboratory:
%
%       http://pscal.ece.gatech.edu/testsys/index.html
%
%   MATPOWER case file data provided by Bruce Wollenberg.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	108	22	0	0	1	1	0	138	1	1.05	0.95;
	2	2	97	20	0	0	1	1	0	138	1	1.05	0.95;
	3	1	180	37	0	0	1	1	0	138	1	1.05	0.95;
	4	1	74	15	0	0	1	1	0	138	1	1.05	0.95;
	5	1	71	14	0	0	1	1	0	138	1	1.05	0.95;
	6	1	136	28	0	-100	2	1	0	138	1	1.05	0.95;
	7	2	125	25	0	0	2	1	0	138	1	1.05	0.95;
	8	1	171	35	0	0	2	1	0	138	1	1.05	0.95;
	9	1	175	36	0	0	1	1	0	138	1	1.05	0.95;
	10	1	195	40	0	0	2	1	0	138	1	1.05	0.95;
	11	1	0	0	0	0	3	1	0	230	1	1.05	0.95;
	12	1	0	0	0	0	3	1	0	230	1	1.05	0.95;
	13	3	265	54	0	0	3	1	0	230	1	1.05	0.95;
	14	2	194	39	0	0	3	1	0	230	1	1.05	0.95;
	15	2	317	64	0	0	4	1	0	230	1	1.05	0.95;
	16	2	100	20	0	0	4	1	0	230	1	1.05	0.95;
	17	1	0	0	0	0	4	1	0	230	1	1.05	0.95;
	18	2	333	68	0	0	4	1	0	230	1	1.05	0.95;
	19	1	181	37	0	0	3	1	0	230	1	1.05	0.95;
	20	1	128	26	0	0	3	1	0	230	1	1.05	0.95;
	21	2	0	0	0	0	4	1	0	230	1	1.05	0.95;
	22	2	0	0	0	0	4	1	0	230	1	1.05	0.95;
	23	2	0	0	0	0	3	1	0	230	1	1.05	0.95;
	24	1	0	0	0	0	4	1	0	230	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf	%	Unit Code
mpc.gen = [
	1	10	0	10	0	1.035	100	1	20	16	0	0	0	0	0	0	0	0	0	0	0;	%	U20
	1	10	0	10	0	1.035	100	1	20	16	0	0	0	0	0	0	0	0	0	0	0;	%	U20
	1	76	0	30	-25	1.035	100	1	76	15.2	0	0	0	0	0	0	0	0	0	0	0;	%	U76
	1	76	0	30	-25	1.035	100	1	76	15.2	0	0	0	0	0	0	0	0	0	0	0;	%	U76
	2	10	0	10	0	1.035	100	1	20	16	0	0	0	0	0	0	0	0	0	0	0;	%	U20 (5)
	2	10	0	10	0	1.035	100	1	20	16	0	0	0	0	0	0	0	0	0	0	0;	%	U20
	2	76	0	30	-25	1.035	100	1	76	15.2	0	0	0	0	0	0	0	0	0	0	0;	%	U76
	2	76	0	30	-25	1.035	100	1	76	15.2	0	0	0	0	0	0	0	0	0	0	0;	%	U76
	7	80	0	60	0	1.025	100	1	100	25	0	0	0	0	0	0	0	0	0	0	0;	%	U100
	7	80	0	60	0	1.025	100	1	100	25	0	0	0	0	0	0	0	0	0	0	0;	%	U100 (10)
	7	80	0	60	0	1.025	100	1	100	25	0	0	0	0	0	0	0	0	0	0	0;	%	U100
	13	95.1	0	80	0	1.02	100	1	197	69	0	0	0	0	0	0	0	0	0	0	0;	%	U197
	13	95.1	0	80	0	1.02	100	1	197	69	0	0	0	0	0	0	0	0	0	0	0;	%	U197
	13	95.1	0	80	0	1.02	100	1	197	69	0	0	0	0	0	0	0	0	0	0	0;	%	U197
	14	0	35.3	200	-50	0.98	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;	%	SynCond (15)
	15	12	0	6	0	1.014	100	1	12	2.4	0	0	0	0	0	0	0	0	0	0	0;	%	U12
	15	12	0	6	0	1.014	100	1	12	2.4	0	0	0	0	0	0	0	0	0	0	0;	%	U12
	15	12	0	6	0	1.014	100	1	12	2.4	0	0	0	0	0	0	0	0	0	0	0;	%	U12
	15	12	0	6	0	1.014	100	1	12	2.4	0	0	0	0	0	0	0	0	0	0	0;	%	U12
	15	12	0	6	0	1.014	100	1	12	2.4	0	0	0	0	0	0	0	0	0	0	0;	%	U12 (20)
	15	155	0	80	-50	1.014	100	1	155	54.3	0	0	0	0	0	0	0	0	0	0	0;	%	U155
	16	155	0	80	-50	1.017	100	1	155	54.3	0	0	0	0	0	0	0	0	0	0	0;	%	U155
	18	400	0	200	-50	1.05	100	1	400	100	0	0	0	0	0	0	0	0	0	0	0;	%	U400
	21	400	0	200	-50	1.05	100	1	400	100	0	0	0	0	0	0	0	0	0	0	0;	%	U400
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50 (25)
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50
	22	50	0	16	-10	1.05	100	1	50	10	0	0	0	0	0	0	0	0	0	0	0;	%	U50 (30)
	23	155	0	80	-50	1.05	100	1	155	54.3	0	0	0	0	0	0	0	0	0	0	0;	%	U155
	23	155	0	80	-50	1.05	100	1	155	54.3	0	0	0	0	0	0	0	0	0	0	0;	%	U155
	23	350	0	150	-25	1.05	100	1	350	140	0	0	0	0	0	0	0	0	0	0	0;	%	U350
];
mpc.gfuIndex = [1 2 5 6 9 10 11 16 17 18 19 20]';% unit index
% add three wind farms
mpc.windfarmIndex = [34,35,36]'; 
addWindFarm = [
    	18	300	0	200	-50	1.05	100	1	300	100	0	0	0	0	0	0	0	0	0	0	0;	%	
    	22	300	0	200	-50	1.05	100	1	300	100	0	0	0	0	0	0	0	0	0	0	0;	%	
    	23	300	0	200	-50	1.05	100	1	300	100	0	0	0	0	0	0	0	0	0	0	0;	%	
    ];
mpc.gen = [mpc.gen; addWindFarm];
% add flag
genFlag = zeros(size(mpc.gen,1),1);
genFlag(mpc.gfuIndex) = 1;
genFlag(mpc.windfarmIndex) = 2;
mpc.gen = [mpc.gen, genFlag];
%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.0026	0.0139	0.4611	175	250	200	0	0	1	-360	360;
	1	3	0.0546	0.2112	0.0572	175	208	220	0	0	1	-360	360;
	1	5	0.0218	0.0845	0.0229	175	208	220	0	0	1	-360	360;
	2	4	0.0328	0.1267	0.0343	175	208	220	0	0	1	-360	360;
	2	6	0.0497	0.192	0.052	175	208	220	0	0	1	-360	360;
	3	9	0.0308	0.119	0.0322	175	208	220	0	0	1	-360	360;
	3	24	0.0023	0.0839	0	400	510	600	1.03	0	1	-360	360;
	4	9	0.0268	0.1037	0.0281	175	208	220	0	0	1	-360	360;
	5	10	0.0228	0.0883	0.0239	175	208	220	0	0	1	-360	360;
	6	10	0.0139	0.0605	2.459	175	193	200	0	0	1	-360	360;
	7	8	0.0159	0.0614	0.0166	175	208	220	0	0	1	-360	360;
	8	9	0.0427	0.1651	0.0447	175	208	220	0	0	1	-360	360;
	8	10	0.0427	0.1651	0.0447	175	208	220	0	0	1	-360	360;
	9	11	0.0023	0.0839	0	400	510	600	1.03	0	1	-360	360;
	9	12	0.0023	0.0839	0	400	510	600	1.03	0	1	-360	360;
	10	11	0.0023	0.0839	0	400	510	600	1.02	0	1	-360	360;
	10	12	0.0023	0.0839	0	400	510	600	1.02	0	1	-360	360;
	11	13	0.0061	0.0476	0.0999	500	600	625	0	0	1	-360	360;
	11	14	0.0054	0.0418	0.0879	500	625	625	0	0	1	-360	360;
	12	13	0.0061	0.0476	0.0999	500	625	625	0	0	1	-360	360;
	12	23	0.0124	0.0966	0.203	500	625	625	0	0	1	-360	360;
	13	23	0.0111	0.0865	0.1818	500	625	625	0	0	1	-360	360;
	14	16	0.005	0.0389	0.0818	500	625	625	0	0	1	-360	360;
	15	16	0.0022	0.0173	0.0364	500	600	625	0	0	1	-360	360;
	15	21	0.0063	0.049	0.103	500	600	625	0	0	1	-360	360;
	15	21	0.0063	0.049	0.103	500	600	625	0	0	1	-360	360;
	15	24	0.0067	0.0519	0.1091	500	600	625	0	0	1	-360	360;
	16	17	0.0033	0.0259	0.0545	500	600	625	0	0	1	-360	360;
	16	19	0.003	0.0231	0.0485	500	600	625	0	0	1	-360	360;
	17	18	0.0018	0.0144	0.0303	500	600	625	0	0	1	-360	360;
	17	22	0.0135	0.1053	0.2212	500	600	625	0	0	1	-360	360;
	18	21	0.0033	0.0259	0.0545	500	600	625	0	0	1	-360	360;
	18	21	0.0033	0.0259	0.0545	500	600	625	0	0	1	-360	360;
	19	20	0.0051	0.0396	0.0833	500	600	625	0	0	1	-360	360;
	19	20	0.0051	0.0396	0.0833	500	600	625	0	0	1	-360	360;
	20	23	0.0028	0.0216	0.0455	500	600	625	0	0	1	-360	360;
	20	23	0.0028	0.0216	0.0455	500	600	625	0	0	1	-360	360;
	21	22	0.0087	0.0678	0.1424	500	600	625	0	0	1	-360	360;
];
%

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [								%	bus	Pmin	Pmax	Qmin	Qmax	Unit Code
	2	1500	0	3	0	        130	    400.6849;	%	1	16	20	0	10	U20 %gfu
	2	1500	0	3	0	        130	    400.6849;	%	1	16	20	0	10	U20 %gfu
	2	1500	0	3	0.014142	16.0811	212.3076;	%	1	15.2	76	-25	30	U76
	2	1500	0	3	0.014142	16.0811	212.3076;	%	1	15.2	76	-25	30	U76
	2	1500	0	3	0	        130	    400.6849;	%	2	16	20	0	10	U20 %gfu %5
	2	1500	0	3	0	        130	    400.6849;	%	2	16	20	0	10	U20 %gfu
	2	1500	0	3	0.014142	16.0811	212.3076;	%	2	15.2	76	-25	30	U76
	2	1500	0	3	0.014142	16.0811	212.3076;	%	2	15.2	76	-25	30	U76
	2	1500	0	3	0.052672	43.6615	781.521;	%	7	25	100	0	60	U100 %gfu
	2	1500	0	3	0.052672	43.6615	781.521;	%	7	25	100	0	60	U100 %gfu %10
	2	1500	0	3	0.052672	43.6615	781.521;	%	7	25	100	0	60	U100 %gfu
	2	1500	0	3	0.00717	    48.5804	832.7575;	%	13	69	197	0	80	U197
	2	1500	0	3	0.00717	    48.5804	832.7575;	%	13	69	197	0	80	U197
	2	1500	0	3	0.00717 	48.5804	832.7575;	%	13	69	197	0	80	U197
	2	1500	0	3	0	        0	0;	%	14					SynCond %15
	2	1500	0	3	0.328412	56.564	86.3852;	%	15	2.4	12	0	6	U12 %gfu
	2	1500	0	3	0.328412	56.564	86.3852;	%	15	2.4	12	0	6	U12 %gfu
	2	1500	0	3	0.328412	56.564	86.3852;	%	15	2.4	12	0	6	U12 %gfu
	2	1500	0	3	0.328412	56.564	86.3852;	%	15	2.4	12	0	6	U12 %gfu
	2	1500	0	3	0.328412	56.564	86.3852;	%	15	2.4	12	0	6	U12 %gfu %20 
	2	1500	0	3	0.008342	12.3883	382.2391;	%	15	54.3	155	-50	80	U155
	2	1500	0	3	0.008342	12.3883	382.2391;	%	16	54.3	155	-50	80	U155
	2	1500	0	3	0.000213	4.4231	395.3749;	%	18	100	400	-50	200	U400
	2	1500	0	3	0.000213	4.4231	395.3749;	%	21	100	400	-50	200	U400
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50 %25
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50
	2	1500	0	3	0	0.001	0.001;	%	22	10	50	-10	16	U50 %30
	2	1500	0	3	0.008342	12.3883	382.2391;	%	23	54.3	155	-50	80	U155
	2	1500	0	3	0.008342	12.3883	382.2391;	%	23	54.3	155	-50	80	U155
	2	1500	0	3	0.004895	11.8495	665.1094;	%	23	140	350	-25	150	U350
];
mpc.gencost([mpc.gfuIndex],[2:7])=0;
addWindFarmCost = [
	2	0	0	0	0	0	0;	%	
	2	0	0	0	0	0	0;	%	
	2	0	0	0	0	0	0;	%	
    ];
mpc.gencost = [mpc.gencost; addWindFarmCost];
% replace the high cost gens
mpc.gencost(12:14,:) = repmat(mpc.gencost(7,:),[3,1]);
%% gas bus data
%1st col: bus index 
%2nd col: bus type: 1 for gas system-in bus (usually a  gas source) 'Q' bus, 2 for gas load constant 'P'(for the actual load is not determined before gas flow）
%bus, 3 for connection bus （doesn't have gas source or load)
%3rd col: Pd(Mm3/day)66.25 %4th col: p0(bar) 
%5th col: pmin %6th col: pmax(bar)
Gbus = [
1   1    0       55.8229     0   77
2   1    0       55.7935     0   77
3   2    3.918   55.6551     30  80
4   3    0       54.1081     0   80 
5   2    0       53.0275     0   77  % gas fired unit as gas load
6   2    4.034   52.2771     30  80
7   2    5.256   52.3726     30  80
8   2    0       59.8520     50  66.2
9   3    0       59.4072     0   66.2 
10  2    6.365   57.5939     30  66.2
11  3    0       56.4185     0   66.2
12  2    2.12    54.5150     0   66.2 %e
13  1    0       53.1879     0   66.2
14  2    0       52.9823     0   66.2
15  2    6.848   51.6530     0   66.2
16  2    15.616  50.00       50  66.2 
17  3    0       55.6233     0   66.2
18  3    0       63.00       0   66.2 
19  2    0.222   35.7445     0   66.2
20  2    1.919   33.8422     25  66.2%e
];

% the gas load of GFU has already been considered in the GEOPF, therefore
% no need to add additional load here
mpc.Gbus=Gbus;
%% gas pipeline data
%1:from bus %2:to bus %3:Cij %4:fmin %5:fmmax(Mm3/day)
Gline = [
1 2   3.012 0 125.08
1 2   3.012 0 125.08
2 3   2.459 0 102.13
2 3   2.459 0 102.13
3 4   1.181 0 49.06
5 6   0.317 0 13.15
6 7   0.386 0 16.01
7 4   0.476 0 19.78
4 14  0.812 0 33.73
8 9   2.694 0 111.88
8 9   0.329 0 13.65
9 10  1.347 0 55.94
9 10  0.164 0 6.83
10 11 1.204 0 50.03
10 11 0.147 0 6.11
11 12 0.929 0 38.60
12 13 0.952 0 39.56
13 14 2.694 0 111.88
14 15 1.905 0 79.11
15 16 1.205 0 50.03
11 17 0.227 0 9.42
17 18 0.80 0 999
18 19 0.041 0 999
19 20 0.167 0 6.93] ; 
% test
% Gline(18,3) = 0.41;
%合并相同节点管道
samepipe=[];
for i=1:size(Gline,1)
    for j=(i+1):size(Gline,1)
        if (Gline(i,1)==Gline(j,1))&&(Gline(i,2)==Gline(j,2))
            Gline(i,3:5)=Gline(i,3:5)+Gline(j,3:5);
            samepipe=[samepipe;j];
        end
    end
end
Gline(samepipe,:)=[]; 
mpc.Gline=Gline;

%% GAS SOURCES
% update the gas load and gas source
%1: located bus %2:Pg %3:Pgmin %4:Pgmax(Mm3/day)
mpc.Gsou=[
1  10.911288    8.87    11.594
2  8.4          0       8.4
5  2.814712     0       4.8
8  22.012       20.34   22.012
13 1.2          0       1.2
14 0.96         0       0.96];
% divide into several small gas sources
% proportion = [5,4,2,1,1,1];
% newGsou = [];
% for i = 1:6
%     addGsou = repmat([Gsou(i,1),Gsou(i,2:4)/proportion(i)],[proportion(i),1]);
%     newGsou = [newGsou; addGsou];
% end
% mpc.Gsou=newGsou;
mpc.gasCompositionForGasSource = [
    91.92 	4.39 	0.53 	0.09 	0.00 	0.76 	2.31 
    86.28 	7.01 	1.21 	0.27 	0.00 	0.50 	4.73 
    91.66 	3.88 	0.46 	0.13 	0.00 	1.54 	2.33 
    92.19 	4.32 	0.43 	0.03 	0.00 	0.76 	2.28 
    97.71 	0.63 	0.07 	0.02 	0.00 	1.12 	0.45 
    94.00 	0.00 	0.00 	0.00 	0.50 	2.50 	2.50  
] / 100;
%% gas-electricity interface %用到GEcon的时候小心，第一列不一定按顺序来
%1st col: gas bus index.  %2nd col: elec bus index
% 用于找已知电力节点的GFU的所在天然气节点
GEcon=[
    4 15
    14 1
    16 2
    18 7
];
mpc.GEcon=GEcon;
%% ptg
% gas bus, electricity bus, effciency, min cap, max cap
mpc.ptg=[
1 18  2 200 0 2 %ptg
5 22  2 200 0 2 %ptg
8 23 2 200 0 2 %ptg
];
%% gas price($/Mm^3)
mpc.Gcost=[
    85000 / 24
    85000 / 24
    85000 / 24
    85000 / 24
    62000 / 24
    62000 / 24 ];

mpc.LCecost=[%NOK挪威克朗/kwh
    0 1/60 1 4 24 %第一行是分段线性函数的时间节点,第二个单位是kw，后面三个单位是kwh
    0 5.6 14.4 10.8 8.8 % Large industry
    0 16.6 70.5 57.1 36.1 % Industry 
    0 18.7 99.6 97.1 56.1 % Commercial
    0 4.2 16.2 11.8 8.6 % Agriculture
    0 0 8.6 8.7 7.4 ];% Residential    
% 转换物理量单位为MWh，货币单位换成成美元
mpc.LCecost(2:end,2) = mpc.LCecost(2:end,2) * 60;
mpc.LCecost(2:end,:) = mpc.LCecost(2:end,:) * 0.1272 * 1000;%阿里汇率20180228

% consumer sectors portions
%bus %large industry %industrial % commercial %agriculture %residential
mpc.consumerSectorPortion = [
    0 0.175  0.2775 0.185  0      0.3625;
    0 0.6529 0.0359 0.0553 0.0218 0.2341;
    0 0.4075 0      0.1175 0      0.4750;
    0 0.2775 0.0925 0.185  0      0.445; 
    0 0      0.1525 0.085  0.37   0.3925;%5
    0 0.175  0.2775 0.185  0      0.3625;
    0 0.6529 0.0359 0.0553 0.0218 0.2341;
    0 0.4075 0      0.1175 0      0.4750;
    0 0.2775 0.0925 0.185  0      0.445; 
    0 0      0.1525 0.085  0.37   0.3925;%10
    0 0.175  0.2775 0.185  0      0.3625;
    0 0.6529 0.0359 0.0553 0.0218 0.2341;
    0 0.4075 0      0.1175 0      0.4750;
    0 0.2775 0.0925 0.185  0      0.445; 
    0 0      0.1525 0.085  0.37   0.3925;%15
    0 0.175  0.2775 0.185  0      0.3625;
    0 0.6529 0.0359 0.0553 0.0218 0.2341;
];

%% ----------gas trasient----------------
% 3 diameter (mm) 4 length (km)
GlineRaw = [   
    1	2	890     4
    1	2	890     4
    2	3	890     6
    2	3	890     6
    3	4	890     26
    5	6	590.1   43
    6	7	590.1	29
    7	4	590.1	19
    4	14	890     55 
    8	9	890     5
    8	9	395     5
    9	10	890     20
    9	10	395.5   20 
    10	11  890     25
    10	11  395.5   25
    11	12  890     42
    12	13  890     40
    13	14  890     5
    14	15  890     10
    15	16  890     25
    11	17  395.5   10.5
    17	18  315.5   26  
    18	19  315.5   98
    19	20  315.5   6
    ];
% measuring unit convert
GlineRaw(:,3) = GlineRaw(:,3) / 1000;
GlineRaw(:,4) = GlineRaw(:,4) * 1000;
%合并相同节点管道
GlineRaw(samepipe,:)=[];
% calculate F
[para] = initializeParameters2();
[~,~,~,~,~,B,~,~,rhon] = unpackPara2(para);
D = GlineRaw(:,3); A = pi*(D/2).^2;
F = sqrt( 4*rhon^2*B^2.*GlineRaw(:,4).*Gline(:,3).^2*(10^6/86400)^2 ./ (D.*A.^2*10^10) );
% F = [19.2846524568081,19.2823964261944,22.9256112527948,22.1074586954129,22.1070928804400,22.0662780083436,22.9256536095011,20.5541972779392,20.5366958945471,20.5294097684344,22.9205329384028,22.9219360363933,22.9333207807464,22.9339748683452,22.9372529249245,21.2722611505256,20.7557112813557,20.6517973350085,20.8138804740731]';
gtd.Gline = [GlineRaw,F];
end

function [para] = initializeParameters2()
% p is the pressure (Pa), Z is the compressibility factor,
% Rspec is the specific gas constant for natural gas (J/kgK), 
% Θis the absolute temperature (K), ρ is the density (kg/m3),
% e the speed of wave propogation in the gas(m/s),
% M is the pipe flow rate (kg/s), 
% x is the distance along the pipeline (m), t is time (s), 
%d is the pipe diameter (m), A is the pipe’s crosssectional area (m2),
%F is the Fanning transmission factor(dimensionless)
% 气体的标准状况的气压温度是有标准规定的
% f=1/F^2, F is Fanning transmission factor, f is friction factor
% 可用稳态公式的Cij推出F: dp2/dx = Q^2/(Cij^2*L)
%%
para.Z = 0.8;
% para.R = 8.314;%J/mol*K
para.R = 8314.51/16;%Rgas
% para.Rgas = para.R/0.016 ;%J/(kg·K),M(CH4) = 16;
para.T = 281.15; %K
para.D = 0.89; % m
% para.D = 0.5; % m
para.A = (para.D/2)^2 * pi; %m2
para.L = 4000;%m pipeLength
para.e = sqrt(para.R*para.Z*para.T);

para.pn = 101325;% Pa, pressure in standard condition(STD) of gas
para.Tn = 298; %K, temerature in STD (in ref<multi...>is set as 288K
para.rhon = para.pn / (para.Z * para.R * para.Tn);
% Csquare = 9.07027/10^10*(10^6/24/3600)^2; % 注意单位换算p和Q的单位, for a representative pipe in ref......
end

function [C_q,R,T,Z_T,Z,B,g,alpha,rhon] = unpackPara2(para)
%UNPACKPARA2 Summary of this function goes here
%   Detailed explanation goes here
%% new
C_q = para.rhon*1000000/86400;
R = para.R;
T = para.T;
Z_T = para.Z; % not sure
Z = para.Z;
B = para.e;
g = 9.98;
alpha = 0;

rhon = para.rhon;
end





