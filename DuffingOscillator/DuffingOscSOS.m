clear 
close all
clc

addpath('cvx\sedumi')  
addpath('cvx\sdpt3')  

%% Variables

mpvar( 'z', [2,1])
mpvar( 'x', [4,1])

% Set variables and decision varaibles in SOS
vars = [z; x]; 
prog = sosprogram(vars);

%% NN Parameters
W1 = [0.9572    -0.8003
    0.4854    -0.1419];
W2 = [0.4218    -0.7922
    0.9157    -0.9595];
W3 = [0.6557    -0.8491];
dim_in = 2; 
dim_hidden = [2,2]; 
dim_out = 1; 
dims = [dim_in, dim_hidden, dim_out];
AF = 'relu'; 

%% Intial region 
% u_min = [-10;-10];
% u_max = [10;10];
u_min = -10^10*ones(2,1);
u_max = 10^10*ones(2,1);

%% IBP
net.weights{1} = W1; net.weights{2} = W2; net.weights{3} = W3;
net.biases{1} = [0;0]; net.biases{2} = [0;0]; net.biases{3} = [0];
net.activation = AF;
net.dims = dims;
[Y_min,Y_max,X_min,X_max,out_min,out_max] = intervalBoundPropagation(u_min,u_max,dim_hidden,net);

%% Constraints
%[eq_constraints, ineq_constraints] = hiddenLayerConstraintsTwoSectors(net,u_min,u_max,z,x);
[eq_constraints, ineq_constraints] = hiddenLayerConstraintsOneSector(net,u_min,u_max,[z(1);z(2)],x);

% Slope constraints
%repeated = 1;
%[eq_rep_constraints,ineq_rep_constraints] = hiddenLayerConstraintsRepeated(net,u_min,u_max,repeated,[z(1);z(2)],x);
v_out = net.weights{end}*x(end - dim_hidden(end) + 1 : end) + net.biases{end};

%% System Dynamics (zeta = 0.5)
dotz1 = z(2);
dotz2 = -z(1) - z(2) - z(1)^3 + v_out;

%% Lyapunov function
orderm = 2;
orderV = 4;
[prog,V] = sospolyvar(prog,monomials([z(1),z(2)],2:orderV),'wscoeff');

phi = 0;
for i = 1:2
    constr = 0.01;%sym(0.01);
    for j = 1:orderV/2
        [prog,eps(i,j)] = sossosvar(prog,1);
        phi = phi+eps(i,j)*vars(i)^(2*j);
        constr = constr-eps(i,j);
    end
    prog = sosineq(prog,-constr);
end

derV = diff(V,z(1))*dotz1 + diff(V,z(2))*dotz2;

%% Psatz Constraints

% dt = 0.1;
% z1plus = z(1) + dt*dotz1;
% z2plus = z(2) + dt*dotz2;
% Vplus = subs(V,[z(1);z(2);],[z1plus;z2plus]);
% derV = Vplus - V;

derVexpr = -derV;

% Ineq constraints
s_NN = cell(size(ineq_constraints,1),1);
for j = 1:size(ineq_constraints,1) 
     if j == 3  || j == 6 || j == 9 || j == 12
     else
     [prog,s_NN{j}] = sossosvar(prog,monomials(vars,0:orderm/2));
     derVexpr = derVexpr - s_NN{j}*ineq_constraints{j,1};
    end
end

% Eq constraints
t_NN = cell(size(eq_constraints,1),1);
for j = 1:size(eq_constraints,1)
    [prog,t_NN{j}] = sospolyvar(prog,monomials(vars,0:orderm));
    derVexpr = derVexpr - t_NN{j}*eq_constraints{j,1};
end

% % Add region as check, (sqrt(2)*u_max(1))^2 - z(1)^2 - z(2)^2;
% ineq_reg1 = (z(1) - u_min(1))*(u_max(1) - z(1));
% ineq_reg2 = (z(2) - u_min(2))*(u_max(2) - z(2));
% [prog,pr1] = sossosvar(prog,monomials(vars,0:orderm/2));
% [prog,pr2] = sossosvar(prog,monomials(vars,0:orderm/2));
%derVexpr = derVexpr - pr1*ineq_reg1 - pr2*ineq_reg2;

%% Solve SOS
prog = sosineq(prog,V - phi);
prog = sosineq(prog,derVexpr);

prog = sossolve(prog);
SOL = sosgetsol(prog,V)

return

%% ROA
syms y_1
xline(u_min(1))
hold on
yline(u_min(2))
xline(u_max(1))
yline(u_max(2))

domain1 = [-15, 15, -15, 15];
for j = 0:50:1000
[C3,h3] = pcontour(SOL,j,domain1,'b',[300, 300]);
end


