function [ fit,imp,arxqs,n4s2 ] = CoupleAnalysis( input,output,inputname,ouputname, startPo, endPo )
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



%plot(time,input,time,output)
%legend(inputname,ouputname)

%Removing Equilibrium Values from the Data
input = input - mean(input);
output = output - mean(output(:,:));
%figure
%plot(time,input,time,output)
%legend(strcat(inputname,':mean'),strcat(ouputname,':mean'))

trainInput = input(startPo:endPo,:);
trainOutput = output(startPo:endPo,:);
testInput = input(startPo:endPo,:);
testOuput = output(startPo:endPo,:);

%{    
if useValidateSet then
    input_size =size(input,1); 
    trainset_size = round(input_size*validateRate);

    trainInput = input(1:trainset_size,:);
    trainOutput = output(1:trainset_size,:);
    testInput = input(trainset_size+1:input_size,:);
    testOuput = output(trainset_size+1:input_size,:);
end_if
%}

Ts = 0.5; % Sampling interval is 0.5 sec
zTrain = iddata(trainOutput,trainInput,Ts);
zTest = iddata(testOuput,testInput,Ts);

% Set train data properties
zTrain.TimeUnit = 'sec';
zTrain.InputName = {inputname};
zTrain.InputUnit = {'Dominance'};
zTrain.OutputName = ouputname;
zTrain.OutputUnit = 'Dominance';

% Set test data properties
zTest.TimeUnit = 'sec';
zTest.InputName = {inputname};
zTest.InputUnit = {'Dominance'};
zTest.OutputName = ouputname;
zTest.OutputUnit = 'Dominance';

%plot(zTrain)
%figure   % Open a new MATLAB Figure window
%plot(zTest) % Plot the validation data

imp = impulseest(zTrain,[-2.5,20],'PW',[],'noncausal');
%spad = spa(zTrain,[],[],[]);
Opt = arxOptions;
%arxqs = arx(zTrain,[4 4 0], Opt);
arxqs = arx(zTrain,[8 8 0], Opt);
Opt2 = n4sidOptions('N4Weight','CVA', 'N4Horizon',[15 29 29]);
n4s2 = n4sid(zTrain, 'best', Opt2);                                 
%compare(zTest,imp,spad,arxqs,n4s2)
[y,fit,x0] = compare(zTest,imp,arxqs,n4s2);
end

