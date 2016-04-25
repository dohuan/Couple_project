% Read the .csv file, output the data for each couple
% 134 couples, each has 8-D feature, total features: 134*8=1072
% longest time stampt: 1232
%%
clear
close all
clc

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
