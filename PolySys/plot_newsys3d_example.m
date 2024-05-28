syms z_1 z_2

SOL_3 =  1.5195*z_1^2 - 1.01*z_1*z_2 + 0.68856*z_2^2;

% MN
domain1 = [-3,3,-3,3];
fmn = fcontour(SOL_3,domain1, 'LineColor', 'k', 'LineWidth', 9);
fmn.LevelList = 4.5;
hold on

set(gca,'LooseInset',get(gca,'TightInset'));
ax2 = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax2,'fontsize',22)
set(gcf,'position',[0,0,(1080+1920)/2,1080])
%legend('NNSOSStability','acZF-1-R','acZF-1')
%lgd = legend({'NNSOSStability','acZF-1-R','acZF-1'},'FontSize',30);
%lgd = legend({'4th Order','acZF-1-R','acZF-1'},'FontSize',30);

xlabel('z_1') 
ylabel('z_2') 

% % PP
% X = [182.0269,    8.6600 ; 8.6600,    6.3167];
% fpp = fcontour([z_1 z_2 ]*X*[z_1;z_2 ],domain1, 'LineColor', 'r', 'LineWidth', 9); %,[-1,1,-1,1]
% %hold on
% fpp.LevelList = 1;

% HY - couldn't find anything suitable
% 
% % PP smallest 
% X = [ 979.0829,   85.1996; 85.1996,   45.4196]
% fpps = fcontour([z_1 z_2 ]*X*[z_1;z_2],domain1, 'LineColor', 'b', 'LineWidth', 9); %,[-1,1,-1,1]
% fpps.LevelList = 1;
% 
% set(gca,'LooseInset',get(gca,'TightInset'));
% ax2 = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax2,'fontsize',22)
% set(gcf,'position',[0,0,(1080+1920)/2,1080])
% %legend('NNSOSStability','acZF-1-R','acZF-1')
% lgd = legend({'NNSOSStability','acZF-1-R','acZF-1'},'FontSize',30);
% %lgd = legend({'4th Order','acZF-1-R','acZF-1'},'FontSize',30);
% 
% xlabel('z_1') 
% ylabel('z_2') 