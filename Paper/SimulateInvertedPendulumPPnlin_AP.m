clear; close all;

load('inv_pend_MPC_bias_free.mat')
T = 100;
for z1i = -3:0.1:3
    for z2i = -10:1:10
        z10 = z1i;
        z20 = z2i;
		z30 = sin(z2i);
		z40 = cos(z2i);
        [tout,zout] = ode15s(@(t,x) pendulum(t,x,W),[0,T],[z10,z20]);
		%plot(tout,zout(:,1));
        plot(zout(:,1),zout(:,2))
		hold on
    end
end

function zdot = pendulum(t,z,W)

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
if T > 1
	T = 1;
elseif T < -1
	T = -1;
end
%T = 0;
%zdot = [z(2); (mass*grav*leng*z(3) - mu*z(2) + T)/(mass*leng^2); z(2)*(z(4)); -z(2)*z(3)];
%T
zdot = [z(2); (mass*grav*leng*sin(z(1)) - mu*z(2) + T)/(mass*leng^2)];

end

function y = ReLU(x)
    y=(x>=0).*x;
end