figure
set(gca,'fontsize',16)

hold on

M = csvread('Ngen300.txt');
plot(mean(M),'LineWidth',2);
xlabel('Time of Evolution');
ylabel('Affinity');
title('The process of affinity maturation for anitbodies');
axis([0 300 0.84 1]);

% % Ngen = 100
% x = 100;
% y = mean(M(:,100)) - 0.002;
% txt = '\uparrow Affinity = 09644';
% text(x,y,txt,'HorizontalAlignment','left');
% 
% % Ngen = 300
% x = 300;
% y = mean(M(:,300)) - 0.003;
% aff = mean(M(:));
% txt = 'Affinity = 0.9721\uparrow ';
% text(x,y,txt,'HorizontalAlignment','right');

hold off