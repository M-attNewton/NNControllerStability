

T = 100;
for z1i = -5:1:5
    for z2i = [-5, 5]
        z10 = z1i;
        z20 = z2i;
        [tout,zout] = ode15s(@(t,x) Duff(t,x),[0,T],[z10,z20]);
		%plot(tout,zout(:,1));
        plot(zout(:,1),zout(:,2),'b','LineWidth',2)
		hold on
    end
end

xlim([-4 4])
ylim([-4 4])
set(gca,'LooseInset',get(gca,'TightInset'));
ax2 = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax2,'fontsize',22)
set(gcf,'position',[0,0,(1080+1920)/2,1080])
%legend('NNSOSStability','acZF-1-R','acZF-1')
%lgd = legend({'NNSOSStability','acZF-1-R','acZF-1'},'FontSize',30);
%lgd = legend({'4th Order','acZF-1-R','acZF-1'},'FontSize',30);

xlabel('z_1') 
ylabel('z_2') 

function zdot = Duff(t,z)
    zeta = 0.5;
    W1 = [0.9572    -0.8003
    0.4854    -0.1419];
    W2 = [0.4218    -0.7922
    0.9157    -0.9595];
    W3 = [0.6557    -0.8491];
    layer1 = ReLU(W1*[z(1);z(2)]);
    layer2 = ReLU(W2*layer1);
    output = W3*layer2;
    %output = 0;

    zdot = [z(2);-z(1) - 2*zeta*z(2) - z(1)^3 + output];
end



function y=ReLU(x)
    y=(x>=0).*x;
end

% return
% T = 100;
% for z1i = -3:0.1:3
%     for z2i = -10:1:10
%         z10 = z1i;
%         z20 = z2i;
% 		z30 = sin(z2i);
% 		z40 = cos(z2i);
%         [tout,zout] = ode15s(@(t,x) pendulum(t,x,W),[0,T],[z10,z20]);
% 		%plot(tout,zout(:,1));
%         plot(zout(:,1),zout(:,2))
% 		hold on
%     end
% end
% 
% function zdot = pendulum(t,z,W)
% 
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
% if T > 1
% 	T = 1;
% elseif T < -1
% 	T = -1;
% end
% %T = 0;
% %zdot = [z(2); (mass*grav*leng*z(3) - mu*z(2) + T)/(mass*leng^2); z(2)*(z(4)); -z(2)*z(3)];
% %T
% zdot = [z(2); (mass*grav*leng*sin(z(1)) - mu*z(2) + T)/(mass*leng^2)];
% 
% end
% 