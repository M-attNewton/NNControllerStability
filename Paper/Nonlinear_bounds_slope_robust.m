clear 
close all
clc

%addpath('cvx\sedumi')  
%addpath('cvx\sdpt3')  

%% Variables

mpvar( 'z', [3,1]) % System variables
mpvar( 'x', [10,1]) % NN variables for [5,5] network
mpvar( 'delta', [1,1]) % Robustness of length

% Set variables and decision varaibles in SOS
vars = [z; x; delta]; 
varsz = z;
varsx = x;
prog = sosprogram(vars);

%% NN Parameters
load('inv_pend_MPC_bias_free.mat')
dim_in = 2; 
dim_hidden = [5,5]; 
dim_out = 1; 
dims = [dim_in, dim_hidden, dim_out];
AF = 'tanh'; 

%% Intial region 
%u_min = -1*ones(dim_in,1);
%u_max = 1*ones(dim_in,1);
%u_min = [-pi;-10];
%u_max = [pi;10];
% u_min = [-0.11;-0.45];
% u_max = [0.11;0.45];

% u_min = [-0.25;-0.75];
% u_max = [0.25;0.75];

% u_min = [-0.4;-1.25];
% u_max = [0.4;1.25];

%u_min = [-0.3;-1.4];
%u_max = [0.3;1.4];

u_min = [-0.1;-0.3];
u_max = [0.1;0.3];


%% IBP on initial region
net.weights{1} = W{1}; net.weights{2} = W{2}; net.weights{3} = W{3};
net.biases{1} = b{1}; net.biases{2} = b{2}; net.biases{3} = b{3};
net.activation = AF;
net.dims = dims;

[Y_min,Y_max,X_min,X_max,out_min,out_max] = intervalBoundPropagation(u_min,u_max,dim_hidden,net);

%% Constraints
%[eq_constraints, ineq_constraints] = hiddenLayerConstraintsTwoSectors(net,u_min,u_max,[z(1);z(2)],x);
[eq_constraints, ineq_constraints] = hiddenLayerConstraintsOneSector(net,u_min,u_max,[z(1);z(2)],x);
% Slope constraints
repeated = 1;
[eq_rep_constraints,ineq_rep_constraints] = hiddenLayerConstraintsRepeated(net,u_min,u_max,repeated,[z(1);z(2)],x);
v_out = net.weights{end}*x(end - dim_hidden(end) + 1 : end) + net.biases{end};
%v_out = 0;
%% System Dynamics
mass = 0.15;
%leng = 0.5;
mu = 0.5;
grav = 9.81;
usat = 1;

dotz1 = z(2); %theta
dotz2 = grav*delta*(z(1) - z(3)) - (mu/mass)*delta^2*z(2) + (1/mass)*delta^2*v_out;
%dotz2 = (mass*grav*delta*(z(1) - z(3)) - mu*z(2) + v_out)/(mass*leng^2); %dot theta
%dotz3 = z(2)*z(4); %sin(theta) z(3) = theta - sin(theta)
%dotz4 = -z(2)*z(3); %cos(theta)

%% Lyapunov function
orderV = 2;
orderLyapV = 2;
% orderx = 0;
% [prog,a1] = sospolyvar(prog,1);
% [prog,a2] = sospolyvar(prog,1);
% [prog,a3] = sospolyvar(prog,1);
% [prog,a4] = sospolyvar(prog,1);
% [prog,a5] = sospolyvar(prog,1);
%a5 = 0;

%V1 = a1*z(2)^2+a2*z(3)^2+a3*z(4)^2+a4*z(4)+a5;
%prog = soseq(prog,a3+a4+a5);
%[prog,epsi1] = sossosvar(prog,1);
%phi = epsi1*(1-z(4));
%[prog,V2] = sospolyvar(prog,monomials([z(1),z(2)],2:orderV),'wscoeff');
%V2 = 0;
%V = V1 + V2;
[prog,V] = sospolyvar(prog,monomials([z(1),z(2)],2:orderLyapV),'wscoeff');


% constr1 = epsi1-0.1;
% for j = 1:orderV/2
% 	%[prog,epsilo1(j)] = sossosvar(prog,1);
% 	%phi = phi+epsilo1(j)*z(1)^(2*j);
% 	epsilo1(j)=0;
% 	constr1 = constr1+epsilo1(j);
% end
% prog = sosineq(prog,constr1);
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

%derV = diff(V,z(2))*dotz2 + diff(V,z(3))*dotz3 + diff(V,z(4))*dotz4;
derV = diff(V,z(1))*dotz1 + diff(V,z(2))*dotz2;
derVexpr = -derV;

% Ineq constraints
s_NN = cell(size(ineq_constraints,1),1);
for j = 1:size(ineq_constraints,1) 
    [prog,s_NN{j}] = sossosvar(prog,monomials(vars,0:orderV/2));
    derVexpr = derVexpr - s_NN{j}*ineq_constraints{j,1};
end

% Eq constraints
t_NN = cell(size(eq_constraints,1),1);
for j = 1:size(eq_constraints,1)
    [prog,t_NN{j}] = sospolyvar(prog,monomials(vars,0:orderV));
    derVexpr = derVexpr - t_NN{j}*eq_constraints{j,1};
end

% Slope constraints
sr_NN = cell(size(ineq_rep_constraints,2),1);
for j = 1:size(ineq_rep_constraints,2) 
    [prog,sr_NN{j}] = sossosvar(prog,monomials(vars,0:orderV/2));
    derVexpr = derVexpr - sr_NN{j}*ineq_rep_constraints{1,j};
end

% % Add sin^2 + cos^2 = 1 constraint
% eq_sys = z(3)^2 + (z(4))^2 - 1;
% [prog,p_sys] = sospolyvar(prog,monomials(vars,0:orderV));
% derVexpr = derVexpr - p_sys*eq_sys;

% Add constraints on z(3)
if abs(u_max(1)) > abs(u_min(1))
    grad_z3 = (u_max(1) - sin(u_max(1)))/abs(u_max(1));
elseif abs(u_max(1)) < abs(u_min(1))
    grad_z3 = (u_min(1) - sin(u_min(1)))/(u_min(1));
else
    grad_z3 = (u_max(1) - sin(u_max(1)))/abs(u_max(1));
end
ineq_sys = z(3)*(grad_z3*z(1) - z(3));
[prog,p_sys] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - p_sys*ineq_sys;

% Add local system constraint
ineq_reg1 = (z(1) - u_min(1))*(u_max(1) - z(1));
ineq_reg2 = (z(2) - u_min(2))*(u_max(2) - z(2));
[prog,pr1] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr2] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - pr1*ineq_reg1 - pr2*ineq_reg2;

% ineq_reg1 = z(1) - u_min(1);
% ineq_reg2 = u_max(1) - z(1);
% ineq_reg3 = z(2) - u_min(2);
% ineq_reg4 = u_max(2) - z(2);
% [prog,pr1] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr2] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr3] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pr4] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - pr1*ineq_reg1 - pr2*ineq_reg2 - pr3*ineq_reg3 - pr4*ineq_reg4;

% ineq_reg5 = 0.2^2 - z(1)^2 - z(2)^2;
% [prog,pr5] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - pr5*ineq_reg5;

% Add input constraints
ineq_input = 1 - v_out^2;
[prog,p_in] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - p_in*ineq_input;

% Robustness constraints
ineq_robust1 = (1/0.2 - delta)*(delta - 1/0.8); %first is lower bound, last is upper bound
[prog,pst1] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - pst1*ineq_robust1;

% % Robustness constraints
% ineq_robust1 = 1/0.51 - delta;
% ineq_robust2 = delta - 1/0.51;
% [prog,pst1] = sossosvar(prog,monomials(vars,0:orderV/2));
% [prog,pst2] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - pst1*ineq_robust1 - pst2*ineq_robust2;

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
