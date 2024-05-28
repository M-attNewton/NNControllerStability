clear; close all;

load('inv_pend_MPC_bias_free.mat')
T = 100;
for z1i = -2:0.1:2
    for z2i = -2:0.5:2
        z10 = z1i;
        z20 = z2i;
        [tout,zout] = ode15s(@(t,x) pendulum(t,x,W,sys),[0,T],[z10,z20]);
		%plot(tout,zout)
        plot(zout(:,1),zout(:,2))
		hold on
    end
end

function zdot = pendulum(t,z,W,sys)

mass = 0.15;
leng = 0.5;
mu = 0.5;
grav = 9.81;
usat = 1;
z = z;
y1 = z(1);
y2 = z(2);
layer1 = tanh(W{1}*[y1;y2]);
layer2 = tanh(W{2}*layer1);
T = W{3}*layer2;
if T > usat
	T = usat;
elseif T < -usat
	T = -usat;
end
%A = [0, 1 ; -grav/leng, -mu/(mass*leng^2) ];
%Adis = eye(2) + 0.02*A;
%T=0;
zdot = [z(2); (mass*grav*leng*z(1) - mu*z(2) + T)/(mass*leng^2)];
%zdot = (sys.A - eye(2))/0.02*[z(1);z(2)] + sys.B*T;
end

function y = ReLU(x)
    y=(x>=0).*x;
end