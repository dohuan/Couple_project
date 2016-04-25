function out = featureExtract()
input_file = 'DataHC.mat'; %Joy stick data
load(input_file);
split_count = 7; % 7
stress_file = 'Globalnew.xls';
couples_count = size(Data,2);% couple data starts at 3rd column , and have 8 elements for each couple
invGain = tf(-1,1,0.5); % inverse static gain
% load stress index info from Globalnew.xls file
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
    if ((size(stress_temp,1)==1)&&(isnan(stress_temp)==0)&&(data_size_temp>300))
        
        savedIndex = savedIndex + 1;
    else
        fprintf(' Discarded\n');
    end
end
end



function [ fit, n4s2 ] = CoupleAnalysis_1( output1,output2,output1name,output2name, startPo, endPo )
% CoupleAnalysis analyze a given couple data  by using 3 models(model_1:imp, model_2: arx,model_3: n4s2)
% It returns the simulated/predicted model output such as fit and output
% parameters.
% Input Arguments:
%  input - input data
%  output - output data
%  inputname - name for input data
%  ouputname - name for output data
%  startPo - start index of input data
%  endPo- end index of input data
% Output Arguments:
%  fit- fit values for 3 models
%  imp - imp model output parameters
%  arxqs- arxqs model output parameters
%  n4s2- n4s2 model output parameters
y = [(output1(startPo:endPo,:)-mean(output1(startPo:endPo,:))) ...
    (output2(startPo:endPo,:)-mean((output2(startPo:endPo,:))))];
Ts = 0.5; % Sampling interval is 0.5 sec
z =  iddata(y,[],Ts);
z.TimeUnit = 'sec';
z.OutputName = {output1name, output2name};
z.OutputUnit = {'Dominance', 'Dominance'};
Opt2 = n4sidOptions('N4Weight','CVA', 'N4Horizon',[15 29 29]);
n4s2 = n4sid(z, 5, Opt2); % 5
end
