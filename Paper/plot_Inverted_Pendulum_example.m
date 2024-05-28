clear; close all;

load('inv_pend_MPC_bias_free.mat')
T = 100;
for z1i = -3:0.1:3
    for z2i = -10:4:10
        z10 = z1i;
        z20 = z2i;
		z30 = sin(z2i);
		z40 = cos(z2i);
        [tout,zout] = ode15s(@(t,x) pendulum(t,x,W),[0,T],[z10,z20]);
		%plot(tout,zout(:,1));
        plot(zout(:,1),zout(:,2),'b','LineWidth',2)
		hold on
    end
end

syms z_1 z_2 

SOL = 428.7713*z_1^4 + 42.0105*z_1^3*z_2 + 24.7307*z_1^2*z_2^2 + 3.3859*z_1*z_2^3 + 0.39349*z_2^4 - 2.215e-11*z_1^3 + 2.764e-13*z_1^2*z_2 - 1.5382e-13*z_1*z_2^2 + 2.0163e-14*z_2^3 + 0.0034796*z_1^2 + 0.00025775*z_1*z_2 + 7.101e-05*z_2^2;

% MN
domain1 = [-0.4,0.4,-1.5,1.5];
fmn = fcontour(SOL,domain1, 'LineColor', 'k', 'LineWidth', 9);
fmn.LevelList = 1;
hold on

% PP
X = [182.0269,    8.6600 ; 8.6600,    6.3167];
fpp = fcontour([z_1 z_2 ]*X*[z_1;z_2 ],domain1, 'LineColor', 'r', 'LineWidth', 9); %,[-1,1,-1,1]
%hold on
fpp.LevelList = 1;

% HY - couldn't find anything suitable

% PP smallest 
X = [ 979.0829,   85.1996; 85.1996,   45.4196]
fpps = fcontour([z_1 z_2 ]*X*[z_1;z_2],domain1, 'LineColor', 'y', 'LineWidth', 9); %,[-1,1,-1,1]
fpps.LevelList = 1;

set(gca,'LooseInset',get(gca,'TightInset'));

ax = gca;
ax.FontSize = 22; 
%ax2 = get(gca,'XTickLabel');
%set(gca,'XTickLabel',ax2,'fontsize',22)
set(gcf,'position',[0,0,(1080+1920)/2,1080])
%legend('NNSOSStability','acZF-1-R','acZF-1')
lgd = legend([fmn,fpp,fpps],{'NNSOSStability','acZF-1-R','acZF-1'},'FontSize',30);
%lgd = legend({'4th Order','acZF-1-R','acZF-1'},'FontSize',30);

xlim([-0.35 0.35])
ylim([-1.5 1.5])

xlabel('z_1') 
ylabel('z_2') 


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