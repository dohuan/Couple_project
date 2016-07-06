function [X,y,couple] = featureExtract(dataType, selectFavorite, percentThres)
input_file = 'DataHC.mat'; %Joy stick data
load(input_file);
split_count = 6; % 6 for Domin, 7 for Warm
stress_file = 'Globalnew.xls';
couples_count = size(Data,2);% couple data starts at 3rd column , and have 8 elements for each couple
invGain = tf(-1,1,0.5); % inverse static gain
% --- load stress index info from Globalnew.xls file
[~, id_strings]=xlsread(stress_file, 'sheet1', 'B2:B141');% extract id
id_num = cellfun(@(x) x(3:5), id_strings, 'UniformOutput', false);
all_data=xlsread(stress_file);
stress_numerics = all_data(:,3);% extract stress indice from third column

results = [];
disp('the number of couples:')
disp(couples_count)
savedIndex = 1;
for k=1:couples_count
    fprintf('processing: %d/%d\n',k,couples_count);
    stress_temp = stress_numerics(strcmp(Data(k).coupleID,id_num));
    data_size_temp = size(Data(k).data(:,1),1);
    if ((size(stress_temp,1)==1)&&(isnan(stress_temp)==0)&&(data_size_temp>300))
        if (strcmp(dataType,'domi')==1)
            h = Data(k).data(:,4); % Dominance info
            w = Data(k).data(:,8); % Dominance info
        else
            h = Data(k).data(:,2); % Warmth info
            w = Data(k).data(:,6); % Warmth info
        end
        h= h(~isnan(h),:);
        w= w(~isnan(w),:);
        resultMatrix = zeros(6,1);
        result_param=struct([]) ;
        count=1;
        
        input_size =size(w,1);
        increase_rate = round(input_size/split_count);
        
        startPo = 1;
        endPo = increase_rate;
        
        for step=1:split_count
            [ fit, result_param{count}] = CoupleAnalysis_1( ...
                w,h,'wife','husband',startPo,endPo );
            resultMatrix(1,count)=split_count;
            resultMatrix(2,count)=startPo;
            resultMatrix(3,count)=endPo;
            resultMatrix(4,count)=1;
            resultMatrix(5:6,count)= fit;
            count=count+1;
            
            startPo = endPo+1;
            if step == (split_count-1)
                endPo = input_size;
            else
                endPo = startPo+increase_rate;
            end
        end
        
        if (strcmp(selectFavorite,'husband')==1)
            maxIndex = find(resultMatrix(6,:)==max(resultMatrix(6,:)));
        elseif (strcmp(selectFavorite,'wife')==1)
            maxIndex = find(resultMatrix(5,:)==max(resultMatrix(5,:)));
        else % favorite mean of husband and wife
            maxixtemp = mean(resultMatrix([5 6],:),1);
            maxIndex = find(maxixtemp==max(maxixtemp));
            if (resultMatrix(5,maxIndex)<percentThres||resultMatrix(6,maxIndex)<percentThres)
                maxIndex = NaN;
                fprintf('fit percentage low: Discarded\n')
            end
        end
        if (isnan(maxIndex)~=1)
            fprintf('Fit percentage: %.1f(wife) %.1f(husband)\n',...
                resultMatrix(5,maxIndex),resultMatrix(6,maxIndex));
            n4s2 = result_param{maxIndex};
            couple(savedIndex).id = Data(k).coupleID;
            couple(savedIndex).stressIndex =...
                stress_numerics(strcmp(couple(savedIndex).id,id_num));
            couple(savedIndex).feature = [];
            couple(savedIndex).f_name = {};
            
            couple(savedIndex) = featureFromModel(couple(savedIndex),n4s2);
            savedIndex = savedIndex + 1;
        end
    else
        fprintf(' Discarded\n');
    end
end
couples_count = size(couple,2); % update the number of couples after discard
for i=1:couples_count
    X(i,:) = couple(i).feature;
    y(i,1) = couple(i).stressIndex;
end
y = zscore(y);
for i=1:size(X,2)
    X(:,i) = zscore(X(:,i));
end
% for i=1:size(X,1)
%     X(i,:) = zscore(X(i,:));
% end
[y,ix] = sort(y);
X = X(ix,:);
couple = couple(ix);
end



function [ fit, n4s2 ] = CoupleAnalysis_1( output1,output2,output1name,output2name, startPo, endPo )
y = [(output1(startPo:endPo,:)-mean(output1(startPo:endPo,:))) ...
    (output2(startPo:endPo,:)-mean((output2(startPo:endPo,:))))];
Ts = 0.5; % Sampling interval is 0.5 sec
z =  iddata(y,[],Ts);
z.TimeUnit = 'sec';
z.OutputName = {output1name, output2name};
z.OutputUnit = {'Dominance', 'Dominance'};
Opt2 = n4sidOptions('N4Weight','CVA', 'N4Horizon',[15 29 29]);
n4s2 = n4sid(z, 5, Opt2); % 5
[~,fit,~] = compare(z,n4s2);
end

function couple = featureFromModel(couple, model)
% --- Create a centerized pole (average of eigenvalues on upper half plane)
[V,A] = eig(model.A);
A = diag(A);
ix = [];
for i=1:length(A)
    if (isreal(A(i))==1 || imag(A(i))>0)
        ix = [ix; i];
    end
end
A = mean(A(ix));
V = mean(V(:,ix),2);
C = model.C*V;
tmp = zeros(size(A,1),1);
n4tmp = ss(A,1,C,[], .5);
[freq_tmp,damp_tmp,~] = damp(n4tmp);
dampfreq_tmp = freq_tmp.*sqrt(1-damp_tmp.^2);

% --- Sort modes by time constant \tau
% [wn,xeta] = damp(model);
% tau = 1./(wn.*xeta);
% [~,ix] = sort(tau);
% [V,A] = eig(model.A);
% A = diag(A);
% A = A(ix);
% A = diag(A);
% V = V(:,ix);
% C = model.C*V;
% tmp = zeros(size(A,1),1);
% [freq_tmp,damp_tmp,~] = damp(model);
% dampfreq_tmp = freq_tmp.*sqrt(1-damp_tmp.^2);


% --- Sort modes by eigenvalues
%[A,I] = sort(diag(A),1,'descend');
% [V,A] = eig(model.A);
% [A,I] = sort(diag(A));
% A = diag(A);
% V = V(:,I);
% C = model.C*V;
% tmp = zeros(size(A,1),1);
% [freq_tmp,damp_tmp,~] = damp(model);
% dampfreq_tmp = freq_tmp.*sqrt(1-damp_tmp.^2);

for j=size(tmp,1)
    B = tmp;
    B(j,:) = 1;
    n4tmp = ss(A,B,C,[], .5);
    
    % --- norm
    couple.feature = ...
        [couple.feature norm(n4tmp,2)];
    couple.f_name = ...
        [couple.f_name ['norm-H2-' num2str(j)]];
    
%     couple.feature = ...
%         [couple.feature norm(n4tmp,inf)];
%     couple.f_name = ...
%         [couple.f_name ['norm-Hinf-' num2str(j)]];
    % --- eigvalue
%     couple.feature = ...
%         [couple.feature eigtmp(j)];
%     couple.f_name = ...
%         [couple.f_name ['eig-' num2str(j)]];
    % --- natural freq
%     couple.feature = ...
%         [couple.feature freq_tmp(j)];
%     couple.f_name = ...
%         [couple.f_name ['naFreq-' num2str(j)]];
    % --- damping ratio
    couple.feature = ...
        [couple.feature damp_tmp(j)];
    couple.f_name = ...
        [couple.f_name ['damp-' num2str(j)]];
    % --- damped freq
%     couple.feature = ...
%         [couple.feature dampfreq_tmp(j)];
%     couple.f_name = ...
%         [couple.f_name ['dampedFreq-' num2str(j)]];
    % --- time constant
%     couple.feature = ...
%         [couple.feature 1/(freq_tmp(j)*damp_tmp(j))];
%         %[couple.feature tautmp(j)];
%     couple.f_name = ...
%         [couple.f_name ['tau-' num2str(j)]];
    % --- overshoot
    ostmp = exp(-pi*damp_tmp(j)/sqrt(1-damp_tmp(j)^2))*100;
    couple.feature = ...
        [couple.feature ostmp];
    couple.f_name = ...
        [couple.f_name ['overS-' num2str(j)]];
    % --- real part
%     couple.feature = ...
%         [couple.feature 1/tautmp(j)];
%     couple.f_name = ...
%         [couple.f_name ['realDist-' num2str(j)]];
    
    %         featurevector_r = ...
    %             [featurevector_r norm(n4tmp,inf)];
    %         fname_r = ...
    %             [fname_r ['norm-Hinf-r-' num2str(j)]];
end

% % --- eigenvalue
% eigtmp = abs(A)*B;
% eigtmp = eigtmp(eigtmp~=0);
% %couple.feature = ...
% %    [couple.feature sort((abs(eig(n4tmp))))'];
% featurevector = ...
%     [featurevector eigtmp];
% fname = ...
%     [fname ['eig-r-' num2str(j)]];
% % --- damping ratio
% [freq_temp,damp_temp,~] = damp(n4tmp);
% freq_temp = freq_temp';
% damp_temp = damp_temp';
% featurevector = ...
%     [featurevector freq_temp];
% for i=1:length(freq_temp)
%     fname = ...
%         [fname ['naFreq-r-' num2str(j) '-' num2str(i)]];
% end
% featurevector = ...
%     [featurevector damp_temp];
% for i=1:length(damp_temp)
%     fname = ...
%         [fname ['damp-r-' num2str(j) '-' num2str(i)]];
% end
% % --- driven frequency \omega_d
% featurevector = ...
%     [featurevector freq_temp.*sqrt(1-damp_temp.^2)];
% for i=1:length(freq_temp)
%     fname = ...
%         [fname ['dampFreq-r-' num2str(j) '-' num2str(i)]];
% end

end


