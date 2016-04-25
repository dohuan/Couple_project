close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))
%load ./coupleHC_n4s2.mat
%load ./coupleHC_arx.mat
load coupleHC_n4s2_modal.mat

nt = size(couple,2);
nf = size(X,2);
%--------- scan thr. split to find best Split for linear and LASSO --------
% count = 1;
% splitSpan=80:100;
% for k=splitSpan
%     split = k;
%     index = randperm(nt);
%     X_train = X(index(1:split),:);
%     y_train = y(index(1:split),:);
% 
%     X_test = X(index(split+1:end),:);
%     y_test = y(index(split+1:end),:);
% 
%     %----------------------- linear regression ----------------------------
% %     B=sdpvar(nf,1);
% %     F=[];
% %     solvesdp(F,norm(y_train-X_train*B,2));
% %     y_guess_train = X_train*value(B);
% %     y_guess_test = X_test*value(B);
% %     error_y_test(count) = sum(sqrt((y_guess_test-y_test).^2));
% %     error_y_train(count) = sum(sqrt((y_guess_train-y_train).^2));
% 
%     %----------------------- LASSO regression -----------------------------
%     [FitObj,y_guess_train,y_guess_test] = lassoRegression(X_train,y_train,...
%                     X_test,y_test);
%     error_y_test(count) = sum(sqrt((y_guess_test-y_test).^2));
%     error_y_train(count) = sum(sqrt((y_guess_train.y_guess_train-y_train).^2));
% 
%     count = count+1;
% end


%------------- try linear and LASSO with a fixed Split --------------------
split = 80;
index = randperm(nt);
% X_train = X(index(1:split),:);
% y_train = y(index(1:split),:);
% 
% X_test = X(index(split+1:end),:);
% y_test = y(index(split+1:end),:);

X_train = X(index,:);
y_train = y(index,:);

X_test = X_train;
y_test = y_train;


% Choose type of regression: Linear or LASSO
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



%%
% -------------- test surviving features --------------
% run = 100;
% featureTrack = zeros(nf,1);
% split = 17;
% for i=1:run
%     index = randperm(nt);
%     X_train = X(index(1:split),:);
%     y_train = y(index(1:split),:);
%     
%     X_test = X(index(split+1:end),:);
%     y_test = y(index(split+1:end),:);
%     
%     [FitObj,~,y_guess_test] = lassoRegression(X_train,y_train,...
%         X_test,y_test);
%     for k=1:nf
%         if (FitObj.B_optimal(k)~=0)
%             featureTrack(k,1)=featureTrack(k,1)+1;
%         end
%     end
% end
% survive_index = find(featureTrack>=20);
% index = randperm(nt);
% X_train = X(index(1:split),survive_index);
% y_train = y(index(1:split),:);
% X_test = X(index(split+1:end),survive_index);
% y_test = y(index(split+1:end),:);
% survive_nf = size(survive_index,1);
% B=sdpvar(survive_nf,1);
% F=[];
% solvesdp(F,norm(y_train-X_train*B,2));
% y_guess_train = X_train*value(B);
% y_guess_test = X_test*value(B);
% [feature_index,feature_name] = xlsread('./featureOutlook.xlsx','sheet1','D1:D54');
% featureDropIndex=[];
% if (featureMode(1)==0)
%     featureDropIndex = [featureDropIndex,49,52];
% end
% if (featureMode(2)==0)
%     featureDropIndex = [featureDropIndex,1:8,13:20,25:40,50,53];
% end
% if (featureMode(3)==0)
%     featureDropIndex = [featureDropIndex,9:12,21:24,41:48,51,54];
% end
% feature_name(featureDropIndex) = [];
% figure(1)
% bar(featureTrack);
% set(gca,'XTick',1:nf,'XTickLabel',feature_name);
% rotateXLabels(gca(),90);
% set(gca,'FontSize',10);
% BOpt = value(B);
%%
%----------------------- plot results -------------------------------------
figure(1)
bar(FitObj.B_optimal);

figure(2)
subplot(2,1,1);
hold on;
plot(y_train,'LineWidth',2);
if (reg_mode == 1)
    plot(y_guess_train.y_guess_train,'r--','LineWidth',2);
else
    plot(y_guess_train,'r--','LineWidth',2);
end
title('train fitting');

subplot(2,1,2);
plot(y_test,'LineWidth',2);
hold on;
plot(y_guess_test,'r--','LineWidth',2);
title('test fitting');
legend('true value','predicted value');



