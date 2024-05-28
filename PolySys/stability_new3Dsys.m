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

% Prune network
 W{1} = W{1}(1:5,:);
 W{2} = W{2}(1:5,1:5);
 W{3} = W{3}(1,1:5);
 b{1} = zeros(size(W{1},1),1);
 b{2} = zeros(size(W{2},1),1);
 b{3} = zeros(size(W{3},1),1);

%% Variables

mpvar( 'z', [3,1]) % System variables
mpvar( 'x', [10,1]) % NN variables for [5,5] network

% Set variables and decision varaibles in SOS
vars = [z; x]; 
varsz = z;
varsx = x;
prog = sosprogram(vars);

%% NN Parameters
dim_in = 3; 
dim_hidden = [5,5]; 
dim_out = 1; 
dims = [dim_in, dim_hidden, dim_out];
AF = 'tanh'; 

%% Intial region 
u_min = -3*ones(3,1);
u_max = 3*ones(3,1);
u_min = [-3;-3;-3];
u_max = [3;3;3];

%% IBP on initial region
net.weights{1} = W{1}; net.weights{2} = W{2}; net.weights{3} = W{3};
b{1} = zeros(size(b{1})); b{2} = zeros(size(b{2})); b{3} = zeros(size(b{3})); 
net.biases{1} = b{1}; net.biases{2} = b{2}; net.biases{3} = b{3};
net.activation = AF;
net.dims = dims;

[Y_min,Y_max,X_min,X_max,out_min,out_max] = intervalBoundPropagation(u_min,u_max,dim_hidden,net);

%% Constraints
%[eq_constraints, ineq_constraints] = hiddenLayerConstraintsTwoSectors(net,u_min,u_max,[z(1);z(2)],x);
[eq_constraints, ineq_constraints] = hiddenLayerConstraintsOneSector(net,u_min,u_max,[z(1);z(2);z(3)],x);
% Slope constraints
repeated = 1;
[eq_rep_constraints,ineq_rep_constraints] = hiddenLayerConstraintsRepeated(net,u_min,u_max,repeated,[z(1);z(2);z(3)],x);
v_out = net.weights{end}*x(end - dim_hidden(end) + 1 : end) + net.biases{end};

%% System Dynamics
dotz1 = -z(1) + z(2) - z(3); 
dotz2 = -z(1)*(z(3) + 1) - z(2);
dotz3 = -z(1) + v_out*100;

%% Lyapunov function
orderV = 2;
orderLyapV = 2;

[prog,Vquad] = sospolyvar(prog,monomials([z(1),z(2),z(3)],2:orderLyapV),'wscoeff');

% [prog,a1] = sospolyvar(prog,1);
% Vhigher = a1*z(1)^4;
Vhigher = 0;
V = Vquad + Vhigher;

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

% Ineq constraints
s_NN = cell(size(ineq_constraints,1),1);
for j = 1:size(ineq_constraints,1) 
    [prog,s_NN{j}] = sossosvar(prog,monomials(vars,0:orderV/2)); %[prog,s_NN{j}] = sossosvar(prog,monomials(vars,0:orderV/2));
    derVexpr = derVexpr - s_NN{j}*ineq_constraints{j,1};
end

% % Eq constraints
% t_NN = cell(size(eq_constraints,1),1);
% for j = 1:size(eq_constraints,1)
%     [prog,t_NN{j}] = sospolyvar(prog,monomials(vars,0:orderV));
%     derVexpr = derVexpr - t_NN{j}*eq_constraints{j,1};
% end

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

% % Add local system constraint
ineq_reg1 = z(1) - u_min(1);
ineq_reg2 = u_max(1) - z(1);
ineq_reg3 = z(2) - u_min(2);
ineq_reg4 = u_max(2) - z(2);
ineq_reg5 = z(3) - u_min(3);
ineq_reg6 = u_max(3) - z(3);
[prog,pr1] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr2] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr3] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr4] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr5] = sossosvar(prog,monomials(vars,0:orderV/2));
[prog,pr6] = sossosvar(prog,monomials(vars,0:orderV/2));
derVexpr = derVexpr - pr1*ineq_reg1 - pr2*ineq_reg2 - pr3*ineq_reg3 - pr4*ineq_reg4 - pr5*ineq_reg5 - pr6*ineq_reg6;

% ineq_reg7 = 0.1^2 - z(1)^2 - z(2)^2 - z(3)^2;
% [prog,pr5] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - pr5*ineq_reg5;

% % Add input constraints
% ineq_input = 1 - v_out^2;
% [prog,p_in] = sossosvar(prog,monomials(vars,0:orderV/2));
% derVexpr = derVexpr - p_in*ineq_input;

%% Solve SOS
prog = sosineq(prog,V - phi);
prog = sosineq(prog,derVexpr);
%prog = sosineq(prog,derVexpr,'sparse');
%prog = sosineq(prog,derVexpr,'sparsemultipartite',{z,x});
options.solver = "mosek";
prog = sossolve(prog);%,options);
SOL = sosgetsol(prog,V)


%% ROA
% orderV = 0;
% mpvar( 'gam', [1,1])
% prog2 = sosprogram(z);
% prog2 = sosdecvar(prog2,gam);
% [prog2,qr1] = sossosvar(prog2,monomials(z,0:orderV/2));
% [prog2,qr2] = sossosvar(prog2,monomials(z,0:orderV/2));
% [prog2,qr3] = sossosvar(prog2,monomials(z,0:orderV/2));
% [prog2,qr4] = sossosvar(prog2,monomials(z,0:orderV/2));
% [prog2,qr5] = sossosvar(prog2,monomials(z,0:orderV/2));
% [prog2,qr6] = sossosvar(prog2,monomials(z,0:orderV/2));
% %SOL = 1.4648*z(1)^2 - 0.97158*z(1)*z(2) + 4.1378*z(1)*z(3) + 0.66135*z(2)^2 - 2.5674*z(2)*z(3) + 6.9403*z(3)^2
% ROA_expr = SOL - gam - qr1*ineq_reg1 - qr2*ineq_reg2 - qr3*ineq_reg3 - qr4*ineq_reg4 - qr5*ineq_reg5 - qr6*ineq_reg6;
% prog2 = sosineq(prog2,ROA_expr);
% prog2 = sossetobj(prog2,-gam);
% prog2 = sossolve(prog2);
% SOL_gam = sosgetsol(prog2,gam)

return

syms y_1
xline(u_min(1))
hold on
yline(u_min(2))
xline(u_max(1))
yline(u_max(2))

SOL_3 = subs(SOL,z(3),0)

domain1 = [-5, 5, -5, 5];
for j = 0:0.5:5
[C3,h3] = pcontour(SOL_3,j,domain1,'b',[300, 300]);
end


syms y_1
xline(u_min(1))
hold on
yline(u_min(3))
xline(u_max(1))
yline(u_max(3))

SOL_2 = subs(SOL,z(2),0)

domain1 = [-5, 5, -5, 5];
for j = 0:0.5:5
[C3,h3] = pcontour(SOL_2,j,domain1,'b',[300, 300]);
end

syms y_1
xline(u_min(2))
hold on
yline(u_min(3))
xline(u_max(2))
yline(u_max(3))

SOL_1 = subs(SOL,z(1),0)

domain1 = [-5, 5, -5, 5];
for j = 0:0.5:5
[C3,h3] = pcontour(SOL_1,j,domain1,'b',[300, 300]);
end


return


