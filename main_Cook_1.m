close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%

load DataHC;
start_ix = 10;
y = [Data(1).data(start_ix:end,2) Data(1).data(start_ix:end,6) ...
                             ones(size(Data(1).data(start_ix:end,2),1),1)];
%load(fullfile(matlabroot,'toolbox','ident','iddemos','data','dcmotordata'));
% t = 1:0.5:500;
% 
% y = [sind(t)' t'];
% y = [y ones(size(y,1),1)];
                         
nt = size(y,1);
                         
A = [NaN NaN;NaN NaN];
B = [NaN;NaN];
C = eye(2);
K = zeros(2,2);
D = zeros(2,1);
x0 = [Data(1).data(start_ix,2); Data(1).data(start_ix,6)];
%x0 = y(1,[1 2])';
ms = modstruc(A,B,C,D,K,x0);
th = ms2th(ms,'d',[],[],0.5);

opt = ssestOptions('InitialState',x0,'Focus','simulation','SearchMethod','lsqnonlin');
th = pem(y,th,opt);

% --- Regenerate the signal
u = ones(nt,1);
e = randn(nt,2);
th = sett(th,0.5);
%y_ = idsim([u e],th,x0);
y_ = idsim(u,th,x0);


figure(1)
subplot(2,1,1)
title('husband')
hold on
plot(y(:,1),'LineWidth',2);
plot(y_(:,1),'r--','LineWidth',2);
hold off

subplot(2,1,2)
title('wife')
hold on
plot(y(:,2),'LineWidth',2);
plot(y_(:,2),'r--','LineWidth',2);
hold off


