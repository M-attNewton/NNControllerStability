%clear; close all;

%load('inv_pend_MPC_bias_free.mat')
% W = 0;
% T = 100;
% for z1i = -1:0.1:1
%     for z2i = -1:0.1:1
%         z10 = z1i;
%         z20 = z2i;
%         [tout,zout] = ode15s(@(t,x) newsys(t,x,W,b),[0,T],[z10,z20]);
% 		%plot(tout,zout)
%         plot(zout(:,1),zout(:,2))
% 		hold on
%     end
% end

 %W = 0;

 W{1} = W{1}(1:5,:);
 W{2} = W{2}(1:5,1:5);
 W{3} = W{3}(1,1:5);
 b{1} = zeros(size(W{1},1),1);
 b{2} = zeros(size(W{2},1),1);
 b{3} = zeros(size(W{3},1),1);

T = 3;

for z1i = -3:1:3 %-1:0.5:1
    for z2i = -3:1:3
        for z3i = -3:1:3
        z10 = z1i;
        z20 = z2i;
        z30 = z3i; 
        z1 = zeros(T/0.1,1); z1(1,1) = z1i;
        z2 = z1; z2(1,1) = z2i;
        z3 = z1; z3(1,1) = z3i;
            for t = 1:1:T/0.1
        %[tout,zout] = ode15s(@(t,x) newsys(t,x,W,b),[0,T],[z10,z20,z30]);
		    %plot(tout,zout)
                layer1 = tanh(W{1}*[z1(t);z2(t);z3(t)] + b{1});
                layer2 = tanh(W{2}*layer1 + b{2});
                U = W{3}*layer2 + b{3};
                z1(t+1) = z1(t) + 0.1*(-z1(t) + z2(t) - z3(t));
                z2(t+1) = z2(t) + 0.1*(-z1(t)*(z3(t) + 1) - z2(t));
                z3(t+1) = z3(t) + 0.1*(-z1(t) + U*100);               
            end
            if abs(z1(end))>10 || abs(z2(end))>10 ||abs(z3(end))>10
                falg = 1;
            end
            plot3(z1,z2,z3)
            hold on
        end
    end
end

		       


% function zdot = newsys(t,z,W,b)
% 
% % y1 = z(1);
% % y2 = z(2);
% % layer1 = tanh(W{1}*[y1;y2] + b{1});
% % layer2 = tanh(W{2}*layer1 + b{2});
% % T = W{3}*layer2 + b{3};
% % T = 0;
% 
% 
% % 1
% %zdot = [z(2); T*z(2)^2 - z(1)];
% 
% % 2
% %zdot = [z(2) - z(1)^3; T];
% 
% % 3
% %zdot = [-z(1)*(0.1 + (z(1) + z(2))^2); (T + z(1))*(0.1 + (z(1) + z(2))^2)];
% 
% y1 = z(1);
% y2 = z(2);
% y3 = z(3);
% layer1 = sigmoid(W{1}*[y1;y2;y3] + b{1});
% layer2 = sigmoid(W{2}*layer1 + b{2});
% T = W{3}*layer2 + b{3};
% %T = 0;
% 
% % 4
% zdot = [-z(1) + z(2) - z(3); -z(1)*(z(3) + 1) - z(2); -z(1) + T];
% 
% % 5
% %zdot = [z(1)^3 - z(2); z(3) ;  T];
% 
% end

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
