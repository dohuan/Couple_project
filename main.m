close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
%addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))

%[X,y,couple] = featureExtract('warm', 'both', 30);
[X,y,couple] = featureExtract('domi', 'wife', []);

nt = size(couple,2);
nf = size(X,2);
kfold = 50;
% --- create fold indexes
index = randperm(nt);
spf = floor(nt/kfold); % sample per fold
for i=1:kfold
    fold(i).X_test = X(index(i+(i-1)*spf:i+i*spf),:);
    fold(i).y_test = X(index(i+(i-1)*spf:i+i*spf),:);
end


%------------- try linear and LASSO with a fixed Split --------------------
%split = 80;
index = randperm(nt);
X_train = X(index,:);
y_train = y(index,:);
X_test = X_train;
y_test = y_train;

% --- Choose type of regression: Linear or LASSO
reg_mode = 1; % 0: linear, 1: Lasso
switch(reg_mode)
    case 0
        %%----------------------- linear regression -----------------------
        B = sdpvar(nf,1);
        F=[];
        solvesdp(F,norm(y_train-X_train*B,2));
        y_guess_train = X_train*value(B);
        y_guess_test = X_test*value(B);
        error_y_test = sum(sqrt((y_guess_test-y_test).^2));
        error_y_train = sum(sqrt((y_guess_train-y_train).^2));
        Bopt = value(B);
    case 1
        %----------------------- lasso regression -------------------------

        [FitObj,y_guess_train,y_guess_test] = lassoRegression(X_train,y_train,...
            X_test,y_test);
end



time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);

%----------------------- plot results -------------------------------------
figure(1)
bar(FitObj.B_optimal);
axis tight
box on
ylabel('weights')
xlabel('features')

[~,ix] = sort(y_train); 
y_fit = y_guess_train.y_guess_train(ix);

figure(2)
plot(y,'LineWidth',2);
hold on
plot(y_fit,'r--','LineWidth',2);
hold off
axis tight
box on
legend('true','fitted')
xlabel('couples')
ylabel('stress index')

fprintf('Covariance matrix: \n')
cov([zscore(y_guess_test) y_test])



