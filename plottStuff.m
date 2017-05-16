clc; close all; clear all; 

W_time = importdata('walltime.txt'); 

x = W_time(2:end,1);

su_w = W_time(2:end,2)./W_time(2:end,3:end); % speedup for all

size = size(su_w); 
width = size(2);

figure
plot(x, su_w,'-*','LineWidth',2)

set(gcf,'color','white')
set(gca,'FontSize',16)
xlabel('Order of matrix','fontsize',16)
xticks([120 360 600 840 1080 1440])

ylabel('Speedup','fontsize',16)

legend('4 p','9 p', ...
 '16 p','25 p','36 p','Location','northwest','Orientation','horizontal');
title('Total speedup') 


x =  W_time(1,3:end); 

figure
plot(x, su_w','-*','LineWidth',2)

set(gcf,'color','white')
set(gca,'FontSize',16)
xlabel('Number of processes','fontsize',16)
xticks(x)
ylim([0 10])
xlim([2 45])

ylabel('Speedup','fontsize',16)

legend('120','240', ...
 '360','480','600','720','840','960','1080','1200','320','1440','Location','northeast','Orientation','vertical');
title('Total speedup') 