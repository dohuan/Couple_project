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

ifplot = 0;

%Removing Equilibrium Values from the Data


%output1 = output1 - mean(output1);
%output2 = output2 - mean(output2);

y = [(output1(startPo:endPo,:)-mean(output1(startPo:endPo,:))) ...
               (output2(startPo:endPo,:)-mean((output2(startPo:endPo,:))))];

% trainInput = output1(startPo:endPo,:);
% trainOutput = output(startPo:endPo,:);
% testInput = output1(startPo:endPo,:);
% testOuput = output(startPo:endPo,:);

Ts = 0.5; % Sampling interval is 0.5 sec
z =  iddata(y,[],Ts);
%z =  iddata(y,ones(size(y,1),1),Ts);
z.TimeUnit = 'sec';
z.OutputName = {output1name, output2name};
%z.Output2Name = {output2name};
z.OutputUnit = {'Dominance', 'Dominance'};

%zTrain = iddata(trainOutput,trainInput,Ts);
%zTest = iddata(testOuput,testInput,Ts);

% Set train data properties
%zTrain.TimeUnit = 'sec';
%zTrain.InputName = {output1name};
%zTrain.InputUnit = {'Dominance'};
% zTrain.OutputName = output2name;
% zTrain.OutputUnit = 'Dominance';

% Set test data properties
% zTest.TimeUnit = 'sec';
% zTest.InputName = {output1name};
% zTest.InputUnit = {'Dominance'};
% zTest.OutputName = output2name;
% zTest.OutputUnit = 'Dominance';

%plot(zTrain)
%figure   % Open a new MATLAB Figure window
%plot(zTest) % Plot the validation data

%imp = impulseest(zTrain,[-2.5,20],'PW',[],'noncausal');
%spad = spa(zTrain,[],[],[]);
%Opt = arxOptions;
%arxqs = arx(zTrain,[4 4 0], Opt);
%arxqs = arx(zTrain,[8 8 0], Opt);

Opt2 = n4sidOptions('N4Weight','CVA', 'N4Horizon',[15 29 29]);
%Opt2 = n4sidOptions('Focus','stability','Display','on');
%n4s2 = n4sid(z, 'best', Opt2);  
n4s2 = n4sid(z, 5, Opt2); % 5 
%compare(zTest,imp,spad,arxqs,n4s2)
%[y,fit,x0] = compare(zTest,imp,arxqs,n4s2);

if (ifplot)
    [y_fit,fit,~] = compare(z,n4s2);
    h = figure(1);
    subplot(2,1,1)
    hold on
    plot(y(:,1));
    plot(y_fit(:,1),'--');
    hold off
    subplot(2,1,2)
    hold on
    plot(y(:,2));
    plot(y_fit(:,2),'--');
    hold off
    legend('true','fitted');
else
    [~,fit,~] = compare(z,n4s2);
end
end

