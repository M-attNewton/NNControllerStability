%clear; close all;

%load('inv_pend_MPC_bias_free.mat')
% % W = 0;
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


 %W = 0;
 W{1} = W{1}(1:5,:);
 W{2} = W{2}(1:5,1:5);
 W{3} = W{3}(1,1:5);
 b{1} = zeros(size(W{1},1),1);
 b{2} = zeros(size(W{2},1),1);
 b{3} = zeros(size(W{3},1),1);
T = 10;
for z1i = -5:2:5
    for z2i = -5:2:5
        for z3i = -5:2:5
        z10 = z1i;
        z20 = z2i;
        z30 = z3i;
        [tout,zout] = ode15s(@(t,x) newsys(t,x,W,b),[0,T],[z10,z20,z30]);
		%plot(tout,zout)
        if any(zout(end,1) > 3) || any(zout(end,2) > 3) || any(zout(end,3) > 3)
            %falg = 1;
            %plot3(zout(:,1),zout(:,2),zout(:,3),'r')
            flag1 = size(zout,1);
            flag2 = flag1; flag3 = flag1;
            blah = 0;
            for j = 1:size(zout,1)
                if abs(zout(j,1)) > 5 && blah == 0
                    flag1 = j;
                    blah = 1;
                end
            end
            blah = 0;
            for j = 1:size(zout,1)
                if abs(zout(j,2)) > 5 && blah == 0
                    flag2 = j;
                    blah = 1;
                end
            end
            blah = 0;
            for j = 1:size(zout,1)
                if abs(zout(j,3)) > 5  && blah == 0
                    flag3 = j;
                    blah = 1;
                end
            end
            
            zout_temp = zout(1:min(flag1,min(flag2,flag3)),:);
            plot3(zout_temp(:,1),zout_temp(:,2),zout_temp(:,3),'r','LineWidth',2)
            %z1_temp = all(abs(zout(end,1)) < 3)
            %axis([-3, 3, -3, 3, -3, 3])
            hold on
        else
            plot3(zout(:,1),zout(:,2),zout(:,3),'b','LineWidth',2)
        end
		hold on
        end
    end
end

syms f(z_1,z_2,z_3)
f(z_1,z_2,z_3) =  1.5195*z_1^2 - 1.01*z_1*z_2 + 4.2912*z_1*z_3 + 0.68856*z_2^2 - 2.6733*z_2*z_3 + 7.21*z_3^2 - 4.5;
interval = [-5, 5, -5, 5, -5, 5,];
fi = fimplicit3(f,interval);%,'EdgeColor','b');
%colormap([0 0 0])
%colormap([0 0 0; 0 0 1;  0 1 0; 1 0 0; 1 1 1]);

%fi.XRange = [0 5];
fi.EdgeColor = 'none';
fi.FaceAlpha = 0.75;
xlim([-4 4])
ylim([-4 4])
zlim([-4,4])
shading interp;
%mycolors = [1 0 0; 1 1 0; 0 0 1];
colormap(spring(50));

set(gca,'LooseInset',get(gca,'TightInset'));
ax2 = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax2,'fontsize',22)
set(gcf,'position',[0,0,(1080+1920)/2,1080])
%legend('NNSOSStability','acZF-1-R','acZF-1')
%lgd = legend({'NNSOSStability','acZF-1-R','acZF-1'},'FontSize',30);
%lgd = legend({'4th Order','acZF-1-R','acZF-1'},'FontSize',30);

xlabel('z_1') 
ylabel('z_2') 
zlabel('z_3') 

% if abs(zout(end,1)) > 1 || abs(zout(end,2)) > 1 || abs(zout(end,3)) > 1
%             %falg = 1;
%             plot3(zout(:,1),zout(:,2),zout(:,3),'r')
%             %axis([-3, 3, -3, 3, -3, 3])
%             hold on
%         else
%             plot3(zout(:,1),zout(:,2),zout(:,3),'g')
%         end


% for z1i = 0.25:0.02:0.27
%     for z2i = 0.08:0.02:0.1
%         for z3i = 0.25:0.02:0.27

function zdot = newsys(t,z,W,b)

% y1 = z(1);
% y2 = z(2);
% layer1 = tanh(W{1}*[y1;y2] + b{1});
% layer2 = tanh(W{2}*layer1 + b{2});
% T = W{3}*layer2 + b{3};
%T = 0;


% 1
%zdot = [z(2); T*z(2)^2 - z(1)];

% 2
%zdot = [z(2) - z(1)^3; T];

% 3
%zdot = [-z(1)*(0.1 + (z(1) + z(2))^2); (T + z(1))*(0.1 + (z(1) + z(2))^2)];

y1 = z(1);
y2 = z(2);
y3 = z(3);

% layer1 = tanh(W{1}*[y1;y2;y3] );
% layer2 = tanh(W{2}*layer1 );
% U = W{3}*layer2 ;

layer1 = tanh(W{1}*[y1;y2;y3] + b{1});
layer2 = tanh(W{2}*layer1 + b{2});
U = W{3}*layer2 + b{3};

%U = 0;
%T
% 4
zdot = [-z(1) + z(2) - z(3); -z(1)*(z(3) + 1) - z(2); -z(1) + U*100];

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
