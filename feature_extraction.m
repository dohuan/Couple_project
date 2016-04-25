close all
clear all
clc
tic
readCSV = 0;
if (readCSV==1)
    %% ----------------------- Read csv file ----------------------------------
    % !!! Remember to change save file name below when change to WC
    fileName = 'Joystick_HC_updated_forJC.csv'; % husband conflict

    fid = fopen(fileName);
    f_str = textscan(fid,'%s',1);
    b = cellstr(f_str{1});
    feature_tag = strsplit(b{1},',');
    feature_tag{1} = [];
    emptyCell = cellfun('isempty',feature_tag);
    feature_tag(emptyCell) = [];
    %time_stamp = csvread(fileName,1,1,[1 1 1232 1]);
    fclose(fid);
    for i=1:134
        if (i==1)
            data_temp = csvread(fileName,1,1,[1 1 1232 1+7*i]);
            Data(i).f_tag = feature_tag(1:1+7*i);
        else
            data_temp = csvread(fileName,1,8*i-7,[1 8*i-7 1232 8*i]);
            % --- i+7(i-1) -> i+7(i)
            Data(i).f_tag = feature_tag(8*i-7:8*i);
        end
        for j=1:1232
            if (nnz(data_temp(j,:))==0)
                cutIndex = j;
                break;
            end
        end
        %Data(i).time_stamp = time_stamp(1:cutIndex-1);
        Data(i).data = data_temp(1:cutIndex-1,:);
        Data(i).coupleID = Data(i).f_tag{1}([2 3 4]);
        fprintf('Reading data... %d%% \n',round(i/134*100));
    end
else
    input_file = 'DataHC.mat'; %Joy stick data
    load(input_file);
end
%--------------------------input setting-----------------------------------
%--- Number of splits on input
split_count = 7;

% -------------------------------------------------------------------------
% Data: + f_tag: feature name tag
%       + data: data
%       + coupleID : 3-digit ID number
% -------------------------------------------------------------------------
%------------------------ stress data file --------------------------------
stress_file = 'Globalnew.xls';  % from "Global_for JC.sav"
%%
%----------------------------starting--------------------------------------
% --- load input data

couples_count = size(Data,2);% couple data starts at 3rd column , and have 8 elements for each couple
invGain = tf(-1,1,0.5); % inverse static gain
% load stress index info from Globalnew.xls file
[id_numerics, id_strings]=xlsread(stress_file, 'sheet1', 'B2:B141');% extract id
id_num = cellfun(@(x) x(3:5), id_strings, 'UniformOutput', false);
all_data=xlsread(stress_file);
stress_numerics = all_data(:,3);% extract stress indice from third column

% all information will be saved to this
results = [];
% feature mode: [1 1 1]; 1: ON, 0: OFF; [imp arx n4s2]
featureMode = [0 0 1];
disp('the number of couples:')
disp(couples_count)
savedIndex = 1;
for k=0:couples_count-1
    fprintf('processing: %d/%d\n',k,couples_count-1)
    % ----------------------- checking stressIndex ------------------------
    %name = txt(1,index_h);
    stress_temp = stress_numerics(strcmp(Data(k+1).coupleID,id_num));
    data_size_temp = size(Data(k+1).data(:,1),1);
    if ((size(stress_temp,1)==1)&&(isnan(stress_temp)==0)&&(data_size_temp>300))
        couple(savedIndex).id = Data(k+1).coupleID;
        couple(savedIndex).stressIndex =...
            stress_numerics(strcmp(couple(savedIndex).id,id_num));
        
        %     index_h = 3+(k*8);    % only take Dominance info
        %     index_w = 3+(k*8)+4;  % only take Dominance info
        %     h = D(1:end, index_h);
        %     w = D(1:end, index_w);
        h = Data(k+1).data(:,4); % Dominance info
        w = Data(k+1).data(:,8); % Dominance info
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
        end
        %-------------- find highest fit split section ------------------------
        % --- find max index for imp
        maxIndex_h_w.imp = find(resultMatrix_h_w(5,:)==max(resultMatrix_h_w(5,:)));
        maxIndex_w_h.imp = find(resultMatrix_w_h(5,:)==max(resultMatrix_w_h(5,:)));
        % --- find max index for arxqs
        maxIndex_h_w.arxqs = find(resultMatrix_h_w(6,:)==max(resultMatrix_h_w(6,:)));
        maxIndex_w_h.arxqs = find(resultMatrix_w_h(6,:)==max(resultMatrix_w_h(6,:)));
        % --- find max index for n4s2
        maxIndex_h_w.n4s2 = find(resultMatrix_h_w(7,:)==max(resultMatrix_h_w(7,:)));
        maxIndex_w_h.n4s2 = find(resultMatrix_w_h(7,:)==max(resultMatrix_w_h(7,:)));
        
        fitTrack(savedIndex).h_w.imp = resultMatrix_h_w(5,maxIndex_h_w.imp);
        fitTrack(savedIndex).w_h.imp = resultMatrix_w_h(5,maxIndex_w_h.imp);
        fitTrack(savedIndex).h_w.arxqs = resultMatrix_h_w(6,maxIndex_h_w.arxqs);
        fitTrack(savedIndex).w_h.arxqs = resultMatrix_w_h(6,maxIndex_w_h.arxqs);
        fitTrack(savedIndex).h_w.n4s2 = resultMatrix_h_w(7,maxIndex_h_w.n4s2);
        fitTrack(savedIndex).w_h.n4s2 = resultMatrix_w_h(7,maxIndex_w_h.n4s2);
        
        tf_temp_imp{savedIndex} = series(series(result_param_h_w(maxIndex_h_w.imp).imp,...
            result_param_w_h(maxIndex_w_h.imp).imp),invGain);
        %[couple{k+1}.feature(1,1),~,couple{k+1}.feature(1,2),~] = ...
        %                                               margin(tf_temp_imp{k+1});
        
        tf_temp_arxqs{savedIndex} = series(series(result_param_h_w(maxIndex_h_w.arxqs).arxqs,...
            result_param_w_h(maxIndex_w_h.arxqs).arxqs),invGain);
        %[couple{k+1}.feature(1,3),~,couple{k+1}.feature(1,4),~] = ...
        %                                               margin(tf_temp_arxqs{k+1});
        
        tf_temp_n4s2{savedIndex} = series(series(result_param_h_w(maxIndex_h_w.n4s2).n4s2,...
            result_param_w_h(maxIndex_w_h.n4s2).n4s2),invGain);
        
        %%
        %-------------- extract features here -----------------------------
        %[couple{k+1}.feature(1,5),~,couple{k+1}.feature(1,6),~] = ...
        %                                            margin(tf_temp_n4s2{k+1});
        %couple{k+1}.feature = eig(tf_temp_n4s2{k+1}.a)';
        couple(savedIndex).feature = [];
        
        % --------------------------- poles -------------------------------
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(pole(tf_temp_arxqs{savedIndex})))'];
        end
        if (featureMode(3)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(pole(tf_temp_n4s2{savedIndex})))'];
        end
        
        % ------------------------- eigenvalue ----------------------------
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(eig(tf_temp_arxqs{savedIndex})))'];
        end
        if (featureMode(3)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(eig(tf_temp_n4s2{savedIndex})))'];
        end
        
        % ------------------------ damping ratio --------------------------
        [freq_temp_arxqs,damp_temp_arxqs,~] = damp(tf_temp_arxqs{savedIndex});
        [freq_temp_n4s2,damp_temp_n4s2,~] = damp(tf_temp_n4s2{savedIndex});
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                freq_temp_arxqs',damp_temp_arxqs'];
        end
        if (featureMode(3)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                freq_temp_n4s2',damp_temp_n4s2'];
        end
        
        %     % - gain margin -
        %     % 1: gain margin, 2: freq_gain, 3: phase margin, 4: freq_phase
        %     [margin_temp_imp(1,1),margin_temp_imp(1,2),...
        %         margin_temp_imp(1,3),margin_temp_imp(1,4)] = ...
        %         margin(tf_temp_imp{k+1});
        %
        %     [margin_temp_arxqs(1,1),margin_temp_arxqs(1,2),...
        %         margin_temp_arxqs(1,3),margin_temp_arxqs(1,4)] = ...
        %         margin(tf_temp_arxqs{k+1});
        %
        %     [margin_temp_n4s2(1,1),margin_temp_n4s2(1,2),...
        %         margin_temp_n4s2(1,3),margin_temp_n4s2(1,4)] = ...
        %         margin(tf_temp_n4s2{k+1});
        %     if (featureMode(1)==1)
        %         couple{k+1}.feature = [couple{k+1}.feature,...
        %                                         margin_temp_imp(1,[1,3])];
        %     end
        %     if (featureMode(2)==1)
        %         couple{k+1}.feature = [couple{k+1}.feature,...
        %                                         margin_temp_arxqs(1,[1,3])];
        %     end
        %     if (featureMode(3)==1)
        %         couple{k+1}.feature = [couple{k+1}.feature,...
        %                                         margin_temp_n4s2(1,[1,3])];
        %     end
        savedIndex = savedIndex + 1;
    else
        fprintf(' Discarded\n');
    end
end

couples_count = size(couple,2); % update the number of couples
for i=1:couples_count
    %X(i,:) = abs(couple{i}.feature);
    X(i,:) = couple(i).feature;
    y(i,1) = couple(i).stressIndex;
end

y = zscore(y);
X = zscore(X);

if (featureMode(2)==1 && featureMode(3)~=1)
    save('coupleHC_arx.mat','X','y','couple');
end
if (featureMode(2)~=1 && featureMode(3)==1)
    save('coupleHC_n4s2.mat','X','y','couple');
end
if (featureMode(2)==1 && featureMode(3)==1)
    save('coupleHC_mix.mat','X','y','couple');
end

% nf = size(X,2);
% meanY = nanmean(y);
% 
% % -------------- assign mean to missing values of X and y -----------------
% for i=1:nf
%     meanX = nanmean(X(:,i));
%     for k=1:couples_count
%         if isinf(X(k,i))||isnan(X(k,i))==1
%             X(k,i) = meanX;
%         end
%     end
% end
% for k=1:couples_count
%     if  isinf(y(k))||isnan(y(k))==1
%         y(k,1) = meanY;
%     end
% end
% y_origin = y;
% X_origin = X;
% 
% for i=1:couples_count
%     X(i,:) = (X_origin(i,:)-mean(X_origin(i,:)))/sqrt(var(X_origin(i,:)));
%     y(i,1) = (y_origin(i,1)-mean(y_origin))/sqrt(var(y_origin));
% end
% y=zscore(y);
% X=zscore(X);
% nf = size(X,2); % number of features

