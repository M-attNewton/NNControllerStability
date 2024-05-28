clear
clc
close all

%% Load NN
load('nn_4_tanh_data.mat')
data = nn_4_tanh_data;
dim_in = data(1);
dim_out = data(2);

NN_data = data(6:end);

W{1} = zeros(20,dim_in); 
b{1} = zeros(20,1);

W{2} = zeros(20,20); 
b{2} = zeros(20,1);

W{3} = zeros(1,20); 
b{3} = zeros(1,1);

k = 1;
for ii = 1:20
    for jj = 1:dim_in
        W{1}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{1}(ii) = NN_data(k);
    k = k + 1;
end

for ii = 1:20
    for jj = 1:20
        W{2}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{2}(ii) = NN_data(k);
    k = k + 1;
end

for ii = 1:1
    for jj = 1:20
        W{3}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{3}(ii) = NN_data(k);
    k = k + 1;
end

%% SOS

%% Variables

mpvar( 'z', [3,1]) % System variables
vars = z;
prog = sosprogram(vars);

%% System Dynamics

dotz1 = -z(1) + z(2) - z(3);
dotz2 = -z(1)*(z(3) + 1) - z(2);
dotz3 = -z(1);


%% Lyapunov function
orderLyapV = 2;
orderV = 2;

[prog,V] = sospolyvar(prog,monomials([z(1),z(2),z(3)],2:orderLyapV),'wscoeff');

phi = 0;
constr1 = -0.1;
for j = 1:orderLyapV/2
	[prog,epsilo1(j)] = sossosvar(prog,1);
	phi = phi+epsilo1(j)*z(1)^(2*j);
	constr1 = constr1+epsilo1(j);
end
prog = sosineq(prog,constr1);
constr2 = -0.1;
for j = 1:orderLyapV/2
	[prog,epsilo2(j)] = sossosvar(prog,1);
	phi = phi+epsilo2(j)*z(2)^(2*j);
	constr2 = constr2+epsilo2(j);
end
prog = sosineq(prog,constr2);
constr3 = -0.1;
for j = 1:orderLyapV/2
	[prog,epsilo3(j)] = sossosvar(prog,1);
	phi = phi+epsilo3(j)*z(3)^(2*j);
	constr3 = constr3+epsilo3(j);
end
prog = sosineq(prog,constr3);

derV = diff(V,z(1))*dotz1 + diff(V,z(2))*dotz2 + diff(V,z(3))*dotz3;
derVexpr = -derV;

ineq_reg5 = 0.5^2 - z(1)^2 - z(2)^2 - z(3)^2;
[prog,pr5] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - pr5*ineq_reg5;


% % Add local system constraint
% ineq_reg1 = z(1) - u_min(1);
% ineq_reg2 = u_max(1) - z(1);
% ineq_reg3 = z(2) - u_min(2);
% ineq_reg4 = u_max(2) - z(2);
% [prog,pr1] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr2] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr3] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr4] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - pr1*ineq_reg1 - pr2*ineq_reg2 - pr3*ineq_reg3 - pr4*ineq_reg4;


%% Solve SOS
prog = sosineq(prog,V - phi);
prog = sosineq(prog,derVexpr);
%prog = sosineq(prog,derVexpr,'sparse');
%prog = sosineq(prog,derVexpr,'sparsemultipartite',{z,x});
options.solver = "mosek";
prog = sossolve(prog);%,options);
SOL = sosgetsol(prog,V)

return

%% ROA

syms y_1
xline(u_min(1))
hold on
yline(u_min(2))
xline(u_max(1))
yline(u_max(2))

domain1 = [-0.5, 0.5, -1.5, 1.5];
for j = 0:0.1:1
[C3,h3] = pcontour(SOL,j,domain1,'b',[300, 300]);
end

%return
figure
syms x1 x2

% PP
X = [182.0269,    8.6600 ; 8.6600,    6.3167];
fpp = fcontour([x1 x2]*X*[x1;x2],[-0.5,0.5,-1.5,1.5]); %,[-1,1,-1,1]
hold on
fpp.LevelList = 1;

% HY - couldn't find anything suitable

% PP smallest 
X = [ 979.0829,   85.1996; 85.1996,   45.4196]
fpps = fcontour([x1 x2]*X*[x1;x2],[-0.5,0.5,-1.5,1.5]); %,[-1,1,-1,1]
fpps.LevelList = 1;

% MN
[C3,h3] = pcontour(SOL,1,domain1,'b',[300, 300]);


%%
syms y_1
xline(u_min(1))
hold on
yline(u_min(2))
%fplot(u_min(2))
xline(u_max(1))
yline(u_max(2))
%fplot(u_max(2))

syms z1 z2
%plot_sol = 1.6219*z2^2 + 0.85271*sin(z1)^2 + 0.86545*cos(z1)^2 - 1.0078*cos(z1);
%plot_sol = 1.6221*z2^2 + 0.85486*sin(z1)^2 + 0.86567*cos(z1)^2 - 1.0104*cos(z1)
%plot_sol = 0.39983*z2^2 + 1.7287*sin(z1)^2 + 0.10653*cos(z1)^2 - 2.0154*cos(z1) + 1.9089;
plot_sol = 19.2432*z1^2 + 0.52142*z1*z2 + 0.28942*z2^2;
domain1 = [-1, 1, -1, 1];
fc = fcontour(plot_sol, domain1)
fc.LevelList = [0:0.1:0.5];
%hold on



%% JUNK

%plot_sol = 1.6151*z2^2 + 0.85224*sin(z1)^2 + 0.78278*cos(z1)^2 - 1.0086*cos(z1);
%plot_sol = 1.6183*z2^2 + 0.8534*sin(z1)^2 + 0.83175*cos(z1)^2 - 1.0094*cos(z1)
%plot_sol =  1.6221*z2^2 + 0.85486*sin(z1)^2 + 0.86567*cos(z1)^2 - 1.0104*cos(z1);
%plot_sol = 1.6217*z2^2 + 0.85281*sin(z1)^2 + 0.86527*cos(z1)^2 - 1.008*cos(z1);
%betaa = 1.75; %0.12
%[C1,h1] = pcontour(a2,1,domain1,'r',[300, 300]);
%hold on
%[C2] = fcontour(plot_sol,domain1);
%fc.LevelList = [betaa];
% [C3,h3] = pcontour(SOLderV,1,domain1,'b',[300, 300]);
%h2.LineColor = mycolor('coolblue');
%h2.LineWidth = 1;


%a5=0;
%prog = soseq(prog,a3+a4+a5);
%phi = epsi1*(1-z(4));
%prog = sosineq(prog,Vexpr);
%dotz3 = z(2)*(z(4)); %sin(theta)
%p1 = 0;
%Vexpr = Vexpr - p1*equal1;

% a = (sqrt(2)*u_max)^2 - z(1)^2 - z(2)^2;
%[prog,p1] = sospolyvar(prog,monomials(varsz,0:orderV/2-1));
% 
%a2 = 1.5^2 - (z(1) - 0.6)^2 - z(2)^2;
%a2 = 1 - (z(1))^2/(0.05)^2 - (z(2)^2)/(0.1)^2;
%[prog,p3] = sossosvar(prog,monomials(varsz,0:orderV/2-1));
%p3 = 0;
%Vexpr = Vexpr - p3*ineq2;
%ineq2 = 0.01^2 - (z(1))^2 - (z(2))^2;
%derVexpr = derVexpr - p4*ineq2;
