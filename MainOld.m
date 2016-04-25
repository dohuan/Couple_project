close all
clear all
clc
%input setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%How many splits on input
split_count =10; %7
%input data file
input_file='Chris_data'; %Joy stick data
%stress data file
stress_file='Globalnew.xls'; % from "Global_for JC.sav"

%starting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load input data
load(input_file);
couples_count = size(txt,2)/8;% couple data starts at 3rd column , and have 8 elements for each couple
time_index =strcmp(txt,'time');
time=D(1:end,time_index);

% load stress index info from Globalnew.xls file
[id_numerics, id_strings]=xlsread(stress_file, 'sheet1', 'B2:B141');% extract id 
id_num = cellfun(@(x) x(3:5), id_strings, 'UniformOutput', false);
all_data=xlsread(stress_file);
stress_numerics = all_data(:,3);% extract stress indice from third column

% all information will be saved to this
results = [];

disp('the number of couples:')
disp(couples_count)

% for all couples 
for k=0:couples_count-1,
%for k=0:1,
    index_h = 3+(k*8);
    index_w = 3+(k*8)+4;
    h = D(1:end, index_h);
    w = D(1:end, index_w);
    h= h(~isnan(h),:);
    w= w(~isnan(w),:);
    
    resultMatrix_w_h = zeros(7,1);
    resultMatrix_h_w = zeros(7,1);
    result_param_w_h=struct([]) ;
    result_param_h_w=struct([]) ;
    count=1;
    for i=1:split_count,

        input_size =size(w,1); 
        current_splits = i;
        increase_rate = round(input_size/current_splits);
        
        startPo = 1;
        endPo = increase_rate;
        
        for step=1:current_splits,

           % analysis for wife input and husband output 
           [ fit_w_h ,imp_w_h,arxqs_w_h,n4s2_w_h ] = CoupleAnalysis( ...
               w,h,'wife','husband',startPo,endPo );
           resultMatrix_w_h(1,count)=i;
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
           resultMatrix_h_w(1,count)=i;
           resultMatrix_h_w(2,count)=startPo;
           resultMatrix_h_w(3,count)=endPo;
           resultMatrix_h_w(4,count)=2;
           resultMatrix_h_w(5:7,count)=cell2mat(fit_h_w);   
           
           result_param_h_w(count).imp = imp_h_w;
           result_param_h_w(count).arxqs = arxqs_h_w;
           result_param_h_w(count).n4s2 = n4s2_h_w;
           
           count=count+1;

            startPo = endPo+1;
            if step == (current_splits-1)
                endPo = input_size;
            else
                endPo = startPo+increase_rate;
            end

        end
    end
    
    % save to the file 
    name = txt(1,index_h);
    %xlswrite(strcat('result_',name{1}(2:4),'_w.xlsx'),resultMatrix_w_h);
    %xlswrite(strcat('result_',name{1}(2:4),'_h.xlsx'),resultMatrix_h_w);
   
    % make a struct for each couple
    field1 = 'id';  value1 = name{1}(2:4);    
    result = struct(field1,value1);
        
    % add a stress index to the struct
    stress=stress_numerics(strcmp(name{1}(2:4),id_num));
    result.stressIndex = stress;
        
    % for w_h
    %best fit among 1,2,3 method
    [max_value_w123,max_index_w123]=max(resultMatrix_w_h(5:7,:),[],2);
    [max_value_w,max_index_w]= max(max(resultMatrix_w_h(5:7,:),[],2));
    result.best_fit_w = MakeBestFitStruct( resultMatrix_w_h, max_value_w,max_index_w,max_index_w123);
    result.best_fit_w.params = result_param_w_h(max_index_w123(max_index_w));
    
    %best fit for each 1,2,3 method
    result.best_fit_w1 = MakeBestFitStruct( resultMatrix_w_h, max_value_w123(1),1,max_index_w123);
    result.best_fit_w1.params = result_param_w_h(max_index_w123(1));
    result.best_fit_w2 = MakeBestFitStruct( resultMatrix_w_h, max_value_w123(2),2,max_index_w123);
    result.best_fit_w2.params = result_param_w_h(max_index_w123(2));
    result.best_fit_w3 = MakeBestFitStruct( resultMatrix_w_h, max_value_w123(3),3,max_index_w123);
    result.best_fit_w3.params = result_param_w_h(max_index_w123(3));
    
    %best fit among 1,2 method
    [max_value_w12,max_index_w12]=max(resultMatrix_w_h(5:6,:),[],2);
    [max_value_w12_r,max_index_w12_r]=max(max(resultMatrix_w_h(5:6,:),[],2));
    result.best_fit_w_12 = MakeBestFitStruct( resultMatrix_w_h, max_value_w12_r,max_index_w12_r,max_index_w12);
    result.best_fit_w_12.params = result_param_w_h(max_index_w12(max_index_w12_r));
    
    %best fit among 1,3 method
    [max_value_w13,max_index_w13]=max(resultMatrix_w_h([5,7],:),[],2);
    [max_value_w13_r,max_index_w13_r]=max(max(resultMatrix_w_h([5,7],:),[],2));
    result.best_fit_w_13 = MakeBestFitStruct( resultMatrix_w_h, max_value_w13_r,max_index_w13_r,max_index_w13);
    result.best_fit_w_13.params = result_param_w_h(max_index_w13(max_index_w13_r));
    
    %best fit among 2,3 method
    [max_value_w23,max_index_w23]=max(resultMatrix_w_h([6,7],:),[],2);
    [max_value_w23_r,max_index_w23_r]=max(max(resultMatrix_w_h([6,7],:),[],2));
    result.best_fit_w_23 = MakeBestFitStruct( resultMatrix_w_h, max_value_w23_r,max_index_w23_r,max_index_w23);
    result.best_fit_w_23.params = result_param_w_h(max_index_w23(max_index_w23_r));
     
    % for h_w
    %best fit among 1,2,3 method
    [max_value_h123,max_index_h123]=max(resultMatrix_h_w(5:7,:),[],2);
    [max_value_h,max_index_h]= max(max(resultMatrix_h_w(5:7,:),[],2));
    result.best_fit_h = MakeBestFitStruct( resultMatrix_h_w, max_value_h,max_index_h,max_index_h123);
    result.best_fit_h.params = result_param_h_w(max_index_h123(max_index_h));
    
    %best fit for each 1,2,3 method
    result.best_fit_h1 = MakeBestFitStruct( resultMatrix_h_w, max_value_h123(1),1,max_index_h123);
    result.best_fit_h1.params = result_param_h_w(max_index_h123(1));
    result.best_fit_h2 = MakeBestFitStruct( resultMatrix_h_w, max_value_h123(2),2,max_index_h123);
    result.best_fit_h2.params = result_param_h_w(max_index_h123(2));
    result.best_fit_h3 = MakeBestFitStruct( resultMatrix_h_w, max_value_h123(3),3,max_index_h123);
    result.best_fit_h3.params = result_param_h_w(max_index_h123(3));
    
    %best fit among 1,2 method
    [max_value_h12,max_index_h12]=max(resultMatrix_h_w(5:6,:),[],2);
    [max_value_h12_r,max_index_h12_r]=max(max(resultMatrix_h_w(5:6,:),[],2));
    result.best_fit_h_12 = MakeBestFitStruct( resultMatrix_h_w, max_value_h12_r,max_index_h12_r,max_index_h12);
    result.best_fit_h_12.params = result_param_h_w(max_index_h12(max_index_h12_r));
    
    %best fit among 1,3 method
    [max_value_h13,max_index_h13]=max(resultMatrix_h_w([5,7],:),[],2);
    [max_value_h13_r,max_index_h13_r]=max(max(resultMatrix_h_w([5,7],:),[],2));
    result.best_fit_h_13 = MakeBestFitStruct( resultMatrix_h_w, max_value_h13_r,max_index_h13_r,max_index_h13);
    result.best_fit_h_13.params = result_param_h_w(max_index_h13(max_index_h13_r));
    
    %best fit among 2,3 method
    [max_value_h23,max_index_h23]=max(resultMatrix_h_w([6,7],:),[],2);
    [max_value_h23_r,max_index_h23_r]=max(max(resultMatrix_h_w([6,7],:),[],2));
    result.best_fit_h_23 = MakeBestFitStruct( resultMatrix_h_w, max_value_h23_r,max_index_h23_r,max_index_h23);
    result.best_fit_h_23.params = result_param_h_w(max_index_h23(max_index_h23_r));
    
    
    % for avg
    resultMatrix_avg = (resultMatrix_h_w(5:7,:)+ resultMatrix_w_h(5:7,:))/2;
    resultMatrix_avg = [resultMatrix_h_w(1:4,:); resultMatrix_avg];
    
     %best fit among 1,2,3 method
    [max_value_avg123,max_index_avg123]=max(resultMatrix_avg(5:7,:),[],2);
    [max_value_avg,max_index_avg]= max(max(resultMatrix_avg(5:7,:),[],2));
    result.best_fit_avg = MakeBestFitStruct( resultMatrix_avg, max_value_avg,max_index_avg,max_index_avg123);
    result.best_fit_avg.w_params = result_param_w_h(max_index_avg123(max_index_avg));
    result.best_fit_avg.h_params = result_param_h_w(max_index_avg123(max_index_avg));
    
    %best fit for each 1,2,3 method
    result.best_fit_avg1 = MakeBestFitStruct( resultMatrix_avg, max_value_avg123(1),1,max_index_avg123);
    result.best_fit_avg1.w_params = result_param_w_h(max_index_avg123(1));
    result.best_fit_avg1.h_params = result_param_h_w(max_index_avg123(1));
    result.best_fit_avg2 = MakeBestFitStruct( resultMatrix_avg, max_value_avg123(2),2,max_index_avg123);
    result.best_fit_avg2.w_params = result_param_w_h(max_index_avg123(2));
    result.best_fit_avg2.h_params = result_param_h_w(max_index_avg123(2));
    result.best_fit_avg3 = MakeBestFitStruct( resultMatrix_avg, max_value_avg123(3),3,max_index_avg123);
    result.best_fit_avg3.w_params = result_param_w_h(max_index_avg123(3));
    result.best_fit_avg3.h_params = result_param_h_w(max_index_avg123(3));
    
    %best fit among 1,2 method
    [max_value_avg12,max_index_avg12]=max(resultMatrix_avg(5:6,:),[],2);
    [max_value_avg12_r,max_index_avg12_r]=max(max(resultMatrix_avg(5:6,:),[],2));
    result.best_fit_avg_12 = MakeBestFitStruct( resultMatrix_avg, max_value_avg12_r,max_index_avg12_r,max_index_avg12);
    result.best_fit_avg_12.w_params = result_param_w_h(max_index_avg12(max_index_avg12_r));
    result.best_fit_avg_12.h_params = result_param_h_w(max_index_avg12(max_index_avg12_r));
    
    %best fit among 1,3 method
    [max_value_avg13,max_index_avg13]=max(resultMatrix_avg([5,7],:),[],2);
    [max_value_avg13_r,max_index_avg13_r]=max(max(resultMatrix_avg([5,7],:),[],2));
    result.best_fit_avg_13 = MakeBestFitStruct( resultMatrix_avg, max_value_avg13_r,max_index_avg13_r,max_index_avg13);
    result.best_fit_avg_13.w_params = result_param_w_h(max_index_avg13(max_index_avg13_r));
    result.best_fit_avg_13.h_params = result_param_h_w(max_index_avg13(max_index_avg13_r));
    
    %best fit among 2,3 method
    [max_value_avg23,max_index_avg23]=max(resultMatrix_avg([6,7],:),[],2);
    [max_value_avg23_r,max_index_avg23_r]=max(max(resultMatrix_avg([6,7],:),[],2));
    result.best_fit_avg_23 = MakeBestFitStruct( resultMatrix_avg, max_value_avg23_r,max_index_avg23_r,max_index_avg23);
    result.best_fit_avg_23.w_params = result_param_w_h(max_index_avg23(max_index_avg23_r));
    result.best_fit_avg_23.h_params = result_param_h_w(max_index_avg23(max_index_avg23_r));
    
  
    c_bf_w= struct2cell(result.best_fit_w);
    c_bf_h=struct2cell(result.best_fit_h);
    c_bf_avg=struct2cell(result.best_fit_avg);
    
    c_bf_w1= struct2cell(result.best_fit_w1);
    c_bf_w2= struct2cell(result.best_fit_w2);
    c_bf_w3= struct2cell(result.best_fit_w3);
    c_bf_w12= struct2cell(result.best_fit_w_12);
    c_bf_w13= struct2cell(result.best_fit_w_13);
    c_bf_w23= struct2cell(result.best_fit_w_23);
    
    c_bf_h1=struct2cell(result.best_fit_h1);
    c_bf_h2=struct2cell(result.best_fit_h2);
    c_bf_h3=struct2cell(result.best_fit_h3);
    c_bf_h12=struct2cell(result.best_fit_h_12);
    c_bf_h13=struct2cell(result.best_fit_h_13);
    c_bf_h23=struct2cell(result.best_fit_h_23);
    
    c_bf_avg1=struct2cell(result.best_fit_avg1);
    c_bf_avg2=struct2cell(result.best_fit_avg2);
    c_bf_avg3=struct2cell(result.best_fit_avg3);
    c_bf_avg12=struct2cell(result.best_fit_avg_12);
    c_bf_avg13=struct2cell(result.best_fit_avg_13);
    c_bf_avg23=struct2cell(result.best_fit_avg_23);
    
    

    BF_IW=['I:wife';c_bf_w(1:6)];
    BF_IH=['I:husband';c_bf_h(1:6)];
    BF_AVG=['I:avg';c_bf_avg(1:6)];
   
    BF_IW1=['I:wife,M:1';c_bf_w1(1:6)];
    BF_IW2=['I:wife,M:2';c_bf_w2(1:6)];
    BF_IW3=['I:wife,M:3';c_bf_w3(1:6)];
    BF_IW12=['I:wife,M:12';c_bf_w12(1:6)];
    BF_IW13=['I:wife,M:13';c_bf_w13(1:6)];
    BF_IW23=['I:wife,M:23';c_bf_w23(1:6)];
    
    BF_IH1=['I:husband,M:1';c_bf_h1(1:6)];
    BF_IH2=['I:husband,M:2';c_bf_h2(1:6)];
    BF_IH3=['I:husband,M:3';c_bf_h3(1:6)];
    BF_IH12=['I:husband,M:12';c_bf_h12(1:6)];
    BF_IH13=['I:husband,M:13';c_bf_h13(1:6)];
    BF_IH23=['I:husband,M:23';c_bf_h23(1:6)];
    
    BF_AVG1=['I:avg,M:1';c_bf_avg1(1:6)];
    BF_AVG2=['I:avg,M:2';c_bf_avg2(1:6)];
    BF_AVG3=['I:avg,M:3';c_bf_avg3(1:6)];
    BF_AVG12=['I:avg,M:12';c_bf_avg12(1:6)];
    BF_AVG13=['I:avg,M:13';c_bf_avg13(1:6)];
    BF_AVG23=['I:avg,M:23';c_bf_avg23(1:6)];
    
    
    summary=[BF_IW,BF_IH,BF_AVG,BF_IW1,BF_IW2,BF_IW3,BF_IW12,BF_IW13,BF_IW23,BF_IH1,BF_IH2,BF_IH3,BF_IH12,BF_IH13,BF_IH23,...
        BF_AVG1,BF_AVG2,BF_AVG3,BF_AVG12,BF_AVG13,BF_AVG23];
    
    stressIndex='NaN';
    if ~isnan(result.stressIndex)
       stressIndex=result.stressIndex;
    end
    % save a summary file for each couple
    summary_filename = strcat('summary_',name{1}(1:6),'_',num2str(stressIndex),'.xlsx');
    xlswrite(summary_filename,summary);
    
   % plot
    PlotFits( result , w,h, time);
    
    
    
    % save all information in Results
    results = [results,result];
    
end

