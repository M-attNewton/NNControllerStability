%clear; close all;

%load('inv_pend_MPC_bias_free.mat')
% % W = 0;
start_nodes = 1; %6
num_nodes = 3; %8
W{1} = W{1}(start_nodes:num_nodes,:);
W{2} = W{2}(start_nodes:num_nodes,start_nodes:num_nodes);
W{3} = W{3}(1,start_nodes:num_nodes);
b{1} = b{1}(start_nodes:num_nodes);
b{2} = b{2}(start_nodes:num_nodes);
b{3} = b{3}(1);
T = 30;
for z1i = -3:1:3
    for z2i = -3:1:3
        z10 = z1i;
        z20 = z2i;
        [tout,zout] = ode15s(@(t,x) newsys(t,x,W,b),[0,T],[z10,z20]);
		%plot(tout,zout)
        plot(zout(:,1),zout(:,2))
		hold on
    end
end

 %W = 0;
%  W{1} = W{1}(1:5,:);
%  W{2} = W{2}(1:5,1:5);
%  W{3} = W{3}(1,1:5);
%  b{1} = zeros(size(W{1},1),1);
%  b{2} = zeros(size(W{2},1),1);
%  b{3} = zeros(size(W{3},1),1);
% T = 100;
% for z1i = -2:1:2
%     for z2i = -2:1:2
%         for z3i = -2:1:2
%         z10 = z1i;
%         z20 = z2i;
%         z30 = z3i;
%         [tout,zout] = ode15s(@(t,x) newsys(t,x,W,b),[0,T],[z10,z20,z30]);
% 		%plot(tout,zout)
%         if abs(zout(end,1))>10 || abs(zout(end,2))>10 ||abs(zout(end,3))>10
%             falg = 1;
%         end
%         plot3(zout(:,1),zout(:,2),zout(:,3))
% 		hold on
%         end
%     end
% end

% for z1i = 0.25:0.02:0.27
%     for z2i = 0.08:0.02:0.1
%         for z3i = 0.25:0.02:0.27

function zdot = newsys(t,z,W,b)

y1 = z(1);
y2 = z(2);
layer1 = ReLU(W{1}*[y1;y2] + b{1});
layer2 = layer1;
for j = 2:100  
    layer2 = ReLU(W{2}*layer2 + b{2});
end
T = W{3}*layer2 + b{3};
%T = -0.01*rand;


% 1
%zdot = [z(2); T*z(2)^2 - z(1)];

% 2
%zdot = [z(2) - z(1)^3; T];

% 3
zdot = [-z(1)*(0.1 + (z(1) + z(2))^2); (T + z(1))*(0.1 + (z(1) + z(2))^2)];

% y1 = z(1);
% y2 = z(2);
% y3 = z(3);

% layer1 = tanh(W{1}*[y1;y2;y3] );
% layer2 = tanh(W{2}*layer1 );
% U = W{3}*layer2 ;

% layer1 = tanh(W{1}*[y1;y2;y3] + b{1});
% layer2 = tanh(W{2}*layer1 + b{2});
% U = W{3}*layer2 + b{3};

%U = 0;
%T
% 4
%zdot = [-z(1) + z(2) - z(3); -z(1)*(z(3) + 1) - z(2); -z(1) + U*100];

% 5
%zdot = [z(1)^3 - z(2); z(3) ;  T];

end

function y = ReLU(x)
    y=(x>=0).*x;
end

function y = sigmoid(x)
    y = 1./(1 + exp(-x));
end

%A = [0, 1 ; -grav/leng, -mu/(mass*leng^2) ];
%Adis = eye(2) + 0.02*A;
%T=0;
%zdot = (sys.A - eye(2))/0.02*[z(1);z(2)] + sys.B*T;
% mass = 0.15;
% leng = 0.5;
% mu = 0.5;
% grav = 9.81;
% usat = 1;
% z = z;
% y1 = z(1);
% y2 = z(2);
% layer1 = tanh(W{1}*[y1;y2]);
% layer2 = tanh(W{2}*layer1);
% T = W{3}*layer2;
% if T > usat
% 	T = usat;
% elseif T < -usat
% 	T = -usat;
% end

%% one stable, two unstable, three unstable, four unstable, five unstable
