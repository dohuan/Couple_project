close all
clear all
clc
%--------------------------input setting-----------------------------------
%How many splits on input
split_count =7;
%input data file
input_file='Chris_data'; %Joy stick data
%stress data file
stress_file='Globalnew.xls';  % from "Global_for JC.sav"

%----------------------------starting--------------------------------------
% load input data
load(input_file);
couples_count = size(txt,2)/8;% couple data starts at 3rd column , and have 8 elements for each couple
time_index =strcmp(txt,'time');
time=D(1:end,time_index);
invGain = tf(-1,1,0.5); % inverse static gain
% load stress index info from Globalnew.xls file
[id_numerics, id_strings]=xlsread(stress_file, 'sheet1', 'B2:B141');% extract id 
id_num = cellfun(@(x) x(3:5), id_strings, 'UniformOutput', false);
all_data=xlsread(stress_file);
stress_numerics = all_data(:,3);% extract stress indice from third column

% all information will be saved to this
results = [];

disp('the number of couples:')
disp(couples_count)

for k=0:couples_count-1
    disp(sprintf('processing: %d/%d',k,couples_count-1))
    index_h = 3+(k*8);    % only take Dominance info
    index_w = 3+(k*8)+4;  % only take Dominance info
    h = D(1:end, index_h); 
    w = D(1:end, index_w);
    h= h(~isnan(h),:);
    w= w(~isnan(w),:);
    
    resultMatrix_w_h = zeros(7,1);
    resultMatrix_h_w = zeros(7,1);
    result_param_w_h=struct([]) ;
    result_param_h_w=struct([]) ;
    count=1;
    
    input_size =size(w,1); 
    %current_splits = i;
    increase_rate = round(input_size/split_count);

    startPo = 1;
    endPo = increase_rate;

    for step=1:split_count

       % analysis for wife input and husband output 
       [ fit_w_h ,imp_w_h,arxqs_w_h,n4s2_w_h ] = CoupleAnalysis( ...
           w,h,'wife','husband',startPo,endPo );
       resultMatrix_w_h(1,count)=split_count;
       resultMatrix_w_h(2,count)=startPo;
       resultMatrix_w_h(3,count)=endPo;
       resultMatrix_w_h(4,count)=1;
       resultMatrix_w_h(5:7,count)=cell2mat(fit_w_h);   

       result_param_w_h(count).imp = imp_w_h;
       result_param_w_h(count).arxqs = arxqs_w_h;
       result_param_w_h(count).n4s2 = n4s2_w_h;

       % analysis for husband input and wife output
       [ fit_h_w ,imp_h_w,arxqs_h_w,n4s2_h_w ] = CoupleAnalysis( ...
           h,w,'husband','wife',startPo,endPo );
       resultMatrix_h_w(1,count)=split_count;
       resultMatrix_h_w(2,count)=startPo;
       resultMatrix_h_w(3,count)=endPo;
       resultMatrix_h_w(4,count)=2;
       resultMatrix_h_w(5:7,count)=cell2mat(fit_h_w);   

       result_param_h_w(count).imp = imp_h_w;
       result_param_h_w(count).arxqs = arxqs_h_w;
       result_param_h_w(count).n4s2 = n4s2_h_w;

       count=count+1;

        startPo = endPo+1;
        if step == (split_count-1)
            endPo = input_size;
        else
            endPo = startPo+increase_rate;
        end
        %tf_temp = series(series(n4s2_w_h,n4s2_h_w),invGain);
        tf_temp = series(series(imp_w_h,imp_h_w),invGain);
        %tf_temp = series(series(arxqs_w_h,arxqs_h_w),invGain);
        [couple{k+1}.feature(2*step-1),~,couple{k+1}.feature(2*step),~] = ...
                                                    margin(tf_temp);
    end
    name = txt(1,index_h);
    couple{k+1}.id = name{1}(2:4);
    couple{k+1}.stressIndex = stress_numerics(strcmp(name{1}(2:4),id_num));
    
end

for i=1:couples_count
   X(i,:) = zscore(couple{i}.feature);
   y(i,1) = couple{i}.stressIndex;
end

nf = size(X,2);
meanY = nanmean(y);
for i=1:nf
   meanX = nanmean(X(:,i));
   for k=1:couples_count
      if isinf(X(k,i))||isnan(X(k,i))==1
          X(k,i) = meanX;
      end
   end
end
for k=1:couples_count
   if  isinf(y(k))||isnan(y(k))==1
       y(k,1) = meanY;
   end
end
y=zscore(y);
% [B, FitInfo] = lasso(X,y);
% index = 57;
% Bopt = B(:,index);
% IntceptOpt = FitInfo.Intercept(index);
% y_guess = X*Bopt + IntceptOpt;
setSplit = 16;
% T=(1:couples_count);
% even = T(mod(T,2)==0);
% odd = T(mod(T,2)~=0);
% 
% X_train = X(even,:);
% X_test = X(odd,:);
% y_train = y(even);
% y_test = y(odd);
X_train = X(1:setSplit-1,:);
X_test = X(setSplit:end,:);
y_train = y(1:setSplit-1,:);
y_test = y(setSplit:end,:);

%[FitObj,y_train_eval,y_guess_test]=lassoRegression(X_train,y_train,X_test,y_test,...
%                                                0.3);
[FitObj,y_train_eval,y_guess_test]=lassoRegression(X_train,y_train,X_test,y_test);

FitObj.mse
figure('name','test')
hold on;
plot(y_guess_test,'bo-');
plot(y_test,'ro--');
hold off;
grid on
figure('name','train');
plot(y_train_eval.y_guess_train,'bo-');
hold on;
plot(y_train_eval.y_train,'ro--');
hold off;
grid on

