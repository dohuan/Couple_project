close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%

load DataHC;
start_ix = 10;

% --- Read stress index from 'Globalnew.xls'
stress_file = 'Globalnew.xls';  % from "Global_for JC.sav"
[id_numerics, id_strings]=xlsread(stress_file, 'sheet1', 'B2:B141');% extract id
id_num = cellfun(@(x) x(3:5), id_strings, 'UniformOutput', false);
all_data=xlsread(stress_file);
stress_numerics = all_data(:,3);% extract stress indice from third column


eigvalues = zeros(2,length(Data));
dist_to_zero = zeros(2,length(Data));
stress_ix = zeros(1,length(Data));
for t=1:length(Data)
    fprintf(['t=' num2str(t) '\n']);
    y = [Data(t).data(start_ix:end,2) Data(t).data(start_ix:end,6) ...
                             ones(size(Data(t).data(start_ix:end,2),1),1)];
    results = fitting(y);
    eigvalues(:,t) = results.eig;
    dist_to_zero(1,t) = norm(eigvalues(1,t));
    dist_to_zero(2,t) = norm(eigvalues(2,t));
    
    h = figure('Position', [200, 200, 1200, 800]);

    subplot(2,1,1)
    title('husband')
    hold on
    plot(y(:,1),'LineWidth',2);
    plot(results.y_fit(:,1),'r--','LineWidth',2);
    hold off
    box on
    
    subplot(2,1,2)
    title('wife')
    hold on
    plot(y(:,2),'LineWidth',2);
    plot(results.y_fit(:,2),'r--','LineWidth',2);
    hold off
    box on
    
    saveas(h,['./save_fig/' num2str(t) '.eps'],'epsc');
    saveas(h,['./save_fig/' num2str(t) '.jpg'],'jpg');
    
    ix = strcmp(id_num,Data(t).coupleID);
    temp = stress_numerics(ix);
    if(isempty(temp)~=1)
        stress_ix(t) = temp;
    else
        stress_ix(t) = NaN;
    end
    
    close(h);
end

[stress_ix,ix] = sort(stress_ix); 
eigvalues = eigvalues(:,ix);

% figure(1)
% hold on
% plot(stress_ix,'k:','LineWidth',2);
% plot(real(eigvalues(1,:)),imag(eigvalues(1,:)),'r*','LineWidth',2);
% plot(real(eigvalues(2,:)),imag(eigvalues(2,:)),'r*','LineWidth',2);
% hold off

figure(1)
hold on
plot(stress_ix,'k:','LineWidth',2);
plot(real(eigvalues(1,:)),'r-','LineWidth',2);
plot(imag(eigvalues(1,:)),'b-','LineWidth',2);
plot(real(eigvalues(2,:)),'r--','LineWidth',2);
plot(imag(eigvalues(2,:)),'b--','LineWidth',2);
hold off

figure(2)
hold on
plot(stress_ix,'k:','LineWidth',2);
plot(dist_to_zero(1,:),'r--','LineWidth',2);
plot(dist_to_zero(2,:),'b-','LineWidth',2);
