close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
%addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))

[X,y,couple] = featureExtract('warm', 'both', 30); % domi 40 warm 30
%[X,y,couple] = featureExtract('warm', 'wife', []);

nt = size(couple,2);
nf = size(X,2);
%kfold = floor(0.3*nt);
kfold = nt;
fprintf('Number of folds: %d\n',kfold);
% --- create fold indexes
index = randperm(nt);
spf = floor(nt/kfold); % sample per fold


% ---- beta-regression
% --- rescale y to be in range (0,1)
a = 0.01;
b = 0.99;
y = a+(b-a)*(y-min(y))/(max(y)-min(y));

R2 = zeros(kfold,1);
vP = zeros(kfold,nf+1);
pre = 0;
for i=1:kfold
    fprintf('fold number: %d/%d\n',i,kfold);
    test_ix = index(i+(i-1)*(spf-1):i+i*(spf-1));
    X_test = X(test_ix,:);
    y_test = y(test_ix,:);
    X_train = X;
    y_train = y;
    X_train(test_ix,:) = [];
    y_train(test_ix,:) = [];
    
	[R2(i),vP(i,:)] = betareg(y_train, X_train);
    etahat = X_test*vP(i,1:end-1)';
    y_pred = exp(etahat) ./ (1+exp(etahat)); 
    pre = pre + mean((y_pred-y_test).^2);
end

fprintf('Mean of R2: %.2f\n',mean(R2));
fprintf('Std of R2: %.2f\n',std(R2));
fprintf('Mean of beta: \n');
mean(vP)
fprintf('PRESS: %.2f\n',pre);




