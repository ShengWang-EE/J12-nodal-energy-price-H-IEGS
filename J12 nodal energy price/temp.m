x = sdpvar(2,1);
y = sdpvar(2,1);
z = sdpvar(2,1);

cons = [ cone([z,x,y]);
    x>=0;
    y >= 0;
    z >= 0;];
objfcn = x+y+z;
options = sdpsettings('solver','gurobi','gurobi.qcpdual',1);
information = optimize(cons,objfcn,options);
dualvar = dual(cons(1));

%%
x = sdpvar(2,1);
Model = [x'*x <= 1];
options = sdpsettings('solver','gurobi','gurobi.qcpdual',1);
optimize(Model,sum(x),options);
dual(Model(1))

%%
Pptg / eta.electrolysis/GCV.hy*3600*24 - QptgMax_hydrogen
%%
windData = xlsread("windspeed_J15.xlsx",1,'E3:G8762')
%%
electricityLoad = xlsread("yearlyElectricityDemandData.csv",1,'D2:D17521');
gasLoadRaw = xlsread("yearlyGasDemandData.csv",1,'H2:H17521');
gasLoad = zeros(3650,1);
for i = 1:365
    gasLoad((i-1)*48+1:i*48) = sum(gasLoadRaw((i-1)*10+1:i*10));
end 
%%
weights = 1:273;
G=graph(mpc.branch(:,1),mpc.branch(:,2),weights);

plot(G,'EdgeLabel',G.Edges.Weight)