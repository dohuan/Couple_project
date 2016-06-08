close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
%addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))

[X,y,couple] = featureExtract('warm', 'both', 30);
%[X,y,couple] = featureExtract('warm', 'wife', []);

nt = size(couple,2);
nf = size(X,2);
%kfold = floor(0.3*nt);
kfold = 40;
fprintf('Number of folds: %d\n',kfold);
% --- create fold indexes
index = randperm(nt);
spf = floor(nt/kfold); % sample per fold

% Btrack = [];
% R2track = [];
% yguessTrack = [];
% ytrueTrack = [];
% minIndextrack =[];

% for i=1:kfold
%     fprintf('fold number: %d/%d\n',i,kfold);
%     test_ix = index(i+(i-1)*(spf-1):i+i*(spf-1));
%     X_test = X(test_ix,:);
%     y_test = y(test_ix,:);
%     X_train = X;
%     y_train = y;
%     X_train(test_ix,:) = [];
%     y_train(test_ix,:) = [];
%     % --- Choose type of regression: Linear or LASSO
%     reg_mode = 1; % 0: linear, 1: Lasso
%     switch(reg_mode)
%         case 0
%             %%----------------------- linear regression -----------------------
%             B = sdpvar(nf,1);
%             F=[];
%             solvesdp(F,norm(y_train-X_train*B,2));
%             y_guess_train = X_train*value(B);
%             y_guess_test = X_test*value(B);
%             error_y_test = sum(sqrt((y_guess_test-y_test).^2));
%             error_y_train = sum(sqrt((y_guess_train-y_train).^2));
%             Bopt = value(B);
%         case 1
%             %----------------------- lasso regression -------------------------
%             
%             [FitObj,y_guess_train,y_guess_test] = ...
%                    lassoRegression(X_train,y_train,X_test,y_test);
%     end
%     Btrack(i,:) = FitObj.B_optimal;
%     %R2track(i,1) = getR2(y_test, y_guess_test);
%     RMSEtrack(i,1) = rmseCal(y_test, y_guess_test);
%     %yguessTrack(i,:) = y_guess_test;
%     %ytrueTrack(i,:) = y_test;
%     %minIndextrack(i,1) = FitObj.minIndex;
% end

[FitObj,y_guess_train,y_guess_test] = ...
                   lassoRegression(X,y,X,y);

R2 = getR2(y, y_guess_test);
MSE = rmseCal(y, y_guess_test);

f_name = couple(1).f_name;
bar(FitObj.B_optimal);
hold on
ix = find(abs(FitObj.B_optimal)>1.5);
B_select = FitObj.B_optimal(ix);
for j=1:length(ix)
    text(ix(j)+2, B_select(j), f_name{ix(j)},'BackgroundColor',[1 1 1]);
end
               
time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);


 % --- Export figures
% fprintf('Mean and Std of RMSE: %.2f %.2f\n',mean(RMSEtrack),std(RMSEtrack));
% f_name = couple(1).f_name;
% h = figure(1);
% meanB = mean(Btrack);
% bar(meanB);
% hold on
% errorbar(meanB,std(Btrack),'x');
% ix = find(abs(meanB)>0.1);
% B_select = meanB(ix);
% for j=1:length(ix)
%     text(ix(j)+2, B_select(j), f_name{ix(j)},'BackgroundColor',[1 1 1]);
% end
% axis tight
% box on
% ylabel('weights')
% xlabel('features')
%saveas(h,[savePath loadFile{i} '_fname_CV.jpg']);

%----------------------- plot results -------------------------------------
% figure(1)
% bar(FitObj.B_optimal);
% axis tight
% box on
% ylabel('weights')
% xlabel('features')
% 
% [~,ix] = sort(y_train);
% y_fit = y_guess_train.y_guess_train(ix);
% 
% figure(2)
% plot(y,'LineWidth',2);
% hold on
% plot(y_fit,'r--','LineWidth',2);
% hold off
% axis tight
% box on
% legend('true','fitted')
% xlabel('couples')
% ylabel('stress index')
% 
% fprintf('Covariance matrix: \n')
% cov([zscore(y_guess_test) y_test])



