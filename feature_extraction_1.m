close all
clear all
clc
set(0,'defaultFigureColor',[1 1 1])
tic
%% Compute H2 and H-inf as feature vector
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
split_count = 9; % 7

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
        
        resultMatrix = zeros(6,1);
        %resultMatrix_h_w = zeros(7,1);
        result_param=struct([]) ;
        %result_param_h_w=struct([]) ;
        count=1;
        
        input_size =size(w,1);
        %current_splits = i;
        increase_rate = round(input_size/split_count);
        
        startPo = 1;
        endPo = increase_rate;
        
        for step=1:split_count
            
            % analysis for wife input and husband output
            [ fit, n4s2] = CoupleAnalysis_1( ...
                w,h,'wife','husband',startPo,endPo );
            
            resultMatrix(1,count)=split_count;
            resultMatrix(2,count)=startPo;
            resultMatrix(3,count)=endPo;
            resultMatrix(4,count)=1;
            %resultMatrix(5:6,count)=cell2mat(fit);
            resultMatrix(5:6,count)= fit;
            
            %result_param_w_h(count).imp = imp_w_h;
            %result_param_w_h(count).arxqs = arxqs_w_h;
            result_param(count).n4s2 = n4s2;
            
            % analysis for husband input and wife output
%             [ fit_h_w ,imp_h_w,arxqs_h_w,n4s2_h_w ] = CoupleAnalysis_1( ...
%                 h,w,'husband','wife',startPo,endPo );
%             resultMatrix_h_w(1,count)=split_count;
%             resultMatrix_h_w(2,count)=startPo;
%             resultMatrix_h_w(3,count)=endPo;
%             resultMatrix_h_w(4,count)=2;
%             resultMatrix_h_w(5:7,count)=cell2mat(fit_h_w);
%             
%             result_param_h_w(count).imp = imp_h_w;
%             result_param_h_w(count).arxqs = arxqs_h_w;
%             result_param_h_w(count).n4s2 = n4s2_h_w;
            
            
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
        %maxIndex_h_w.imp = find(resultMatrix_h_w(5,:)==max(resultMatrix_h_w(5,:)));
        %maxIndex_w_h.imp = find(resultMatrix(5,:)==max(resultMatrix(5,:)));
        % --- find max index for arxqs
        %maxIndex_h_w.arxqs = find(resultMatrix_h_w(6,:)==max(resultMatrix_h_w(6,:)));
        %maxIndex_w_h.arxqs = find(resultMatrix(6,:)==max(resultMatrix(6,:)));
        % --- find max index for n4s2
        %maxIndex_h_w.n4s2 = find(resultMatrix_h_w(7,:)==max(resultMatrix_h_w(7,:)));
        %maxIndex_w_h.n4s2 = find(resultMatrix(7,:)==max(resultMatrix(7,:)));
        
        maxixtemp = mean(resultMatrix([5 6],:),2);
        maxIndex = find(maxixtemp==max(maxixtemp));
        
        fprintf('Fit percentage: %.1f %.1f\n',resultMatrix(5,maxIndex),resultMatrix(6,maxIndex));
%         fitTrack(savedIndex).h_w.imp = resultMatrix_h_w(5,maxIndex_h_w.imp);
%         fitTrack(savedIndex).w_h.imp = resultMatrix(5,maxIndex.imp);
%         fitTrack(savedIndex).h_w.arxqs = resultMatrix_h_w(6,maxIndex_h_w.arxqs);
%         fitTrack(savedIndex).w_h.arxqs = resultMatrix(6,maxIndex.arxqs);
%         fitTrack(savedIndex).h_w.n4s2 = resultMatrix_h_w(7,maxIndex_h_w.n4s2);
        fitTrack(savedIndex,:) = resultMatrix([5 6],maxIndex);
        
        
        
%         tf_temp_imp{savedIndex} = series(series(result_param_h_w(maxIndex_h_w.imp).imp,...
%             result_param(maxIndex.imp).imp),invGain);
        %[couple{k+1}.feature(1,1),~,couple{k+1}.feature(1,2),~] = ...
        %                                               margin(tf_temp_imp{k+1});
        
%         tf_temp_arxqs{savedIndex} = series(series(result_param_h_w(maxIndex_h_w.arxqs).arxqs,...
%             result_param(maxIndex.arxqs).arxqs),invGain);
        %[couple{k+1}.feature(1,3),~,couple{k+1}.feature(1,4),~] = ...
        %                                               margin(tf_temp_arxqs{k+1});
        
        %tf_temp_n4s2{savedIndex} = series(series(result_param(maxIndex).n4s2,...
        %    result_param(maxIndex).n4s2),invGain);
        
        %n4s2_d = result_param(maxIndex).n4s2;
        %n4s2 = d2c(n4s2_d,'zoh');
        
        n4s2 = result_param(maxIndex).n4s2;
        
        %-------------- extract features here -----------------------------
        %[couple{k+1}.feature(1,5),~,couple{k+1}.feature(1,6),~] = ...
        %                                            margin(tf_temp_n4s2{k+1});
        %couple{k+1}.feature = eig(tf_temp_n4s2{k+1}.a)';
        couple(savedIndex).feature = [];
        %K_track(k,:) = reshape(n4s2.K,[],1);
        % --------------------------- poles -------------------------------
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(pole(tf_temp_arxqs{savedIndex})))'];
        end
        if (featureMode(3)==1)
            %couple(savedIndex).feature = [couple(savedIndex).feature,...
            %    sort(abs(pole(n4s2)))'];
        end
        
        % ------------------------- eigenvalue ----------------------------
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                sort(abs(eig(tf_temp_arxqs{savedIndex})))'];
        end
        if (featureMode(3)==1)
            %couple(savedIndex).feature = [couple(savedIndex).feature,...
            %    sort(abs(eig(n4s2)))'];
            
            %couple(savedIndex).feature = [couple(savedIndex).feature,...
            %    eig(n4s2)'];
            
            tmp = eig(n4s2)';
            if size(tmp~=6)
                tmp(6) = NaN;
            end
            %[tmp,~] = pzmap(n4s2);
            %couple(savedIndex).feature = [couple(savedIndex).feature,...
            %    tmp'];
        end
        
        % ------------------------ damping ratio --------------------------
        %[freq_temp_arxqs,damp_temp_arxqs,~] = damp(tf_temp_arxqs{savedIndex});
        [freq_temp_n4s2,damp_temp_n4s2,~] = damp(n4s2);
        if (featureMode(2)==1)
            couple(savedIndex).feature = [couple(savedIndex).feature,...
                freq_temp_arxqs',damp_temp_arxqs'];
        end
        if (featureMode(3)==1)
            %couple(savedIndex).feature = [couple(savedIndex).feature,...
            %    freq_temp_n4s2',damp_temp_n4s2'];
        end
        
        % --------------------- H2 and H-inf norm -------------------------
%         A = n4s2.C*n4s2.A\n4s2.C;
%         %B = ones(size(A,1),1);
%         %D = n4s2.D;
%         C{1} = [1 0];
%         C{2} = [0 1];
%         C{3} = eye(2);
%         B{1} = [1; 0];
%         B{2} = [0; 1];
%         B{3} = [1; 1];
%         for i=1:length(C)
%             for j=1:length(B)
%                 n4tmp = ss(A,B{j},C{i},[], .5);
%                 couple(savedIndex).feature = ...
%                                 [couple(savedIndex).feature norm(n4tmp,2)];
%                 couple(savedIndex).feature = ...
%                               [couple(savedIndex).feature norm(n4tmp,inf)];
%             end
%         end

        % ----------- extract features with modal decomposition -----------
        [V,A] = eig(n4s2.A);
        %A = n4s2.C*n4s2.A\n4s2.C;
        %B = ones(size(A,1),1);
        %D = n4s2.D;
        C = n4s2.C*V;
        tmp = zeros(size(A,1),1);
        couple(savedIndex).f_name = {};
        for j=1:size(tmp,1)
            B = tmp;
            B(j,:) = 1;
            n4tmp = ss(A,B,C,[], .5);
            % --- H2 and Hinf
            couple(savedIndex).feature = ...
                            [couple(savedIndex).feature norm(n4tmp,2)];
            for i=1:length(norm(n4tmp,2))
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['norm_H2_' num2str(i)]];
            end
            couple(savedIndex).feature = ...
                          [couple(savedIndex).feature norm(n4tmp,inf)];
            for i=1:length(norm(n4tmp,inf))
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['norm_Hinf_' num2str(i)]];
            end
            % --- eigenvalue
            couple(savedIndex).feature = ...
                           [couple(savedIndex).feature (abs(eig(n4tmp)))'];
            eig_track(j,:) = (abs(eig(n4tmp)))';
            for i=1:length(abs(eig(n4tmp)))
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['eig_' num2str(i)]];
            end
            % --- damping ratio
            [freq_temp,damp_temp,~] = damp(n4s2);
            freq_temp = freq_temp';
            damp_temp = damp_temp';
            couple(savedIndex).feature = ...
                        [couple(savedIndex).feature freq_temp];
            for i=1:length(freq_temp)
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['naFreq_' num2str(i)]];
            end
            couple(savedIndex).feature = ...
                        [couple(savedIndex).feature damp_temp];
            for i=1:length(damp_temp)
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['damp_' num2str(i)]];
            end
            % --- driven frequency \omega_d
            couple(savedIndex).feature = ...
              [couple(savedIndex).feature freq_temp.*sqrt(1-damp_temp.^2)];
            
            for i=1:length(freq_temp)
                couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['drivenFreq_' num2str(i)]];
            end
        end
        couple(savedIndex).feature = ...
                        [couple(savedIndex).feature (mean(eig_track,2))'];
        for i=1:length(mean(eig_track,2))
            couple(savedIndex).f_name = ...
                    [couple(savedIndex).f_name ['meanEig_' num2str(i)]];
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

% --- plot K_track
% hold on
% for i=1:size(K_track,2)
%     plot(K_track(:,i));
% end
% hold off

couples_count = size(couple,2); % update the number of couples
for i=1:couples_count
    %X(i,:) = abs(couple{i}.feature);
    X(i,:) = couple(i).feature;
    y(i,1) = couple(i).stressIndex;
end

y = zscore(y);
for i=1:size(X,2)
    X(:,i) = zscore(X(:,i));
end
%X = zscore(X);

% --- plot feature vs target to see if there is any trend 
[y,ix] = sort(y); 
X = X(ix,:);

% --- Compute mean and variance
% upix = find(y>.8);
% stat(1).mean = nanmean([reshape(real(X(upix,:)),[],1) reshape(imag(X(upix,:)),[],1)]);
% stat(1).var = nancov([reshape(real(X(upix,:)),[],1) reshape(imag(X(upix,:)),[],1)]);
% 
% downix = find(y<.2);
% stat(2).mean = nanmean([reshape(real(X(downix,:)),[],1) reshape(imag(X(downix,:)),[],1)]);
% stat(2).var = nancov([reshape(real(X(downix,:)),[],1) reshape(imag(X(downix,:)),[],1)]);
% 
% delix = find(y<.8 & y>.2);
% X(delix,:) = [];
% y(delix)   = [];
% 
% % --- plot root with colors
% colmap = colormap(summer);
% figure(1)
% hold on
% for i=1:size(y,1)
%     coltmp = colmap(round(i/size(y,1)*size(colmap,1)),:);
%     plot(real(X(i,:)),imag(X(i,:)),'o','MarkerFaceColor',coltmp);
% end
% %viscircles([0 0],1);
% error_ellipse('C',stat(1).var,'mu',stat(1).mean','style','k-','conf',.95);
% error_ellipse('C',stat(2).var,'mu',stat(2).mean','style','k--','conf',.95);
% hold off

X_mu = mean(X,2);

subplot(3,1,1)
plot(y)
box on
ylabel('stress index')
axis tight
subplot(3,1,2)
hold on
for i=1:size(X,2)
    plot(X(:,i),'--');
end
hold off
box on
ylabel('H2/Hinf')
axis tight
subplot(3,1,3)
plot(X_mu);
box on
ylabel('features mean')
axis tight

figure(2)
plot(y)
hold on
for i=size(X,2)-5:size(X,2)
    plot(X(:,i),'--');
end
hold off

if (featureMode(2)==1 && featureMode(3)~=1)
    save('coupleHC_arx_single_sys.mat','X','y','couple');
end
if (featureMode(2)~=1 && featureMode(3)==1)
    save('coupleHC_n4s2_single_sys.mat','X','y','couple');
end
if (featureMode(2)==1 && featureMode(3)==1)
    save('coupleHC_mix_single_sys.mat','X','y','couple');
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

