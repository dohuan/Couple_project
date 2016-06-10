%function [bparams,bcost,y_angles,exp_angles,torque,timeout]=fitting_main()
%use either 'odeCody' for Cody's model or 'odefile4' for full state feedback
%update x0 accordingly to be either 9 states with 1st state is 'y' or 8 states with 1st state is 'y'
%update the 'pupper' and 'plower' acoordingly
%update var_ind accordingly (1:8 or 1:7)
clear
clc
tic

addpath(pwd)
prompt = 'What is the model?[1:Cody, 2:Chen02, 3:Peng96,4:Cody2(delay in the loop)]';
modelind = input(prompt);
models={'Cody','Chen02','Peng96','Cody2'};
task={'_nrpt_4_deg','_nfept','_EOnps','_ECnps','_neft','_nfft'};
visit={'_v1','_v2'};
ID={'neck01','neck02','neck03','neck04','neck05',...
    'neck06','neck07','neck08','neck09','neck10'};

%% Optimization
w=1;
options=optimset('maxfunevals',12000,'tolx',1e-5,'tolfun',1e-5,...
    'MaxIter',10000,'Display','iter','PlotFcns',@optimplotresnorm,...%@optimplotfval,...
    'Algorithm','trust-region-reflective');
%'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian'
%'Algorithm','trust-region-reflective'

for j=1:length(ID)
cd(ID{j})

for k=1:1%task index{'_nrpt_4_deg','_nfept','_EOnps','_ECnps','_neft','_nfft'}
    cd1=strcat(ID{j},task{k});
    cd(cd1)
    for z=1:length(visit)
    cd2=strcat(cd1,visit{z});
    cd(cd2);
    s=load(strcat(cd2,'data'),'input','output');%,'startpoints')
    Input=s.input*pi/180;%rad
    Output=s.output*pi/180;%rad
    time=(0:1/60:29.99)';
    if z==1&&j==1
        switch models{modelind}
            case {'Cody','Cody2'}
                %Set the values of the known parameters.  These values are taken from
                %literature, Cody's paper/dissertation
                K=0.605;
                b=0.3;
                I=0.016;
                tc=0.1;
                %Set upper and lower bounds on the fitted parameters.
                pupper=[0.5   10   10   1];
                plower=[0.1  0.05   0.05   0.003];
                const_param=[K,b,I,tc];
                var_ind=1:length(plower); %[delay,Kp,Ki,Kd,k,b,I,tc]
                const_ind=length(plower)+1:length(plower)+length(const_param);
            case 'Chen02'
                K=1;
                pupper=[3   10   0.5 ];%[B, Kvis, delay]
                plower=[0.1  0.05   0.1];
                const_param=[K];
                var_ind=1:length(plower); %[B, Kvis, delay, K]
                const_ind=length(plower)+1:length(plower)+length(const_param);
            case 'Peng96'
                AllInd=1:14;
                var_ind=[1:10,12,13];
                const_ind=sort(setdiff(AllInd,var_ind));
                load('meanParamsPeng','mean_bparamsAllPeng');
                I=0.0148; B=0.1; K=2.077; tc=0.1;%from Cody
                %[Kvis,Kvcr,Kccr,delay,T1a,Tcns1,Tc,Tcns2,Tms1,Tms2,I,B,K,tc]
                pupper=[500 30  0.2  0.4 0.2  1    7  60 1    1    I  2    3  tc];
                plower=[50 0.01 0.01 0.1 0.01 0.04 0.1 2 0.01 0.01 I 0.05 0.5 tc];
                const_param=[I,tc];
                pupper=pupper(var_ind);
                plower=plower(var_ind);
        end
        offset=mean([pupper;plower]);
        gain=pupper-offset;
        pupper2=(pupper-offset)./(gain);
        plower2=(plower-offset)./(gain);
        %% startpoints
%"Optimal Learning" by Warren Powell, rule od thumb is to perform a Latin
%hypercube design with a 2d plus 2 points
        montenum=2*length(plower)+2;
        switch models{modelind}
            case {'Cody','Cody2','Chen02','Peng96'}
        startpoints=find_starts(montenum,plower2,pupper2,var_ind,const_param,const_ind,offset,gain,models{modelind},time);
%         [~,~,~,A]=odeCody(pupper2,var_ind,const_param,const_ind,offset,gain,[],[],[],[],'buildsys',[],[]);
            case ''
                load('startpoints12Peng96','startpoints');
        end
    %% Initialize matrices
    ntrials=size(Input,2);
    %ntrials is usually 3 trials/subject
    bparams1=zeros(montenum,length(pupper));
    bcost1=zeros(montenum,1);
    exitflag1=zeros(montenum,1);
    output1=cell(1,montenum);
    end
    exp_data=Output;
%     x0=zeros(length(A),ntrials);
%     for i=1:ntrials
%         x0(:,i)=[exp_data{i}(1);zeros(8,1)];
%     end
    
    if isempty(gcp('nocreate'))%check if a parallel pool is running or not
        if montenum>16
            parpool(ceil(montenum/2))
        else
        parpool(montenum)
        end
    end
    switch models{modelind}
        case 'Cody'
    parfor i=1:montenum
        %[fun,y,flag,A,B,C]=odeCody(var_param,var_ind,const_param,const_ind,offset,gain,input,time,exp_data,x0,testflag,w,ntrials)
        %[x,resnorm,residual,exitflag,output]=lsqnonlin(fun,x0,lb,ub,options)
        [bparams1(i,:),bcost1(i),~,exitflag1(i),output1{i}] = ...
            lsqnonlin(@(p)odeCody(p,var_ind,const_param,const_ind,offset,gain,Input,time,exp_data,[],'none',w,ntrials),...
            startpoints(i,:),plower2,pupper2,options);
    end
        case 'Chen02'
    parfor i=1:montenum
        %[fun,y,flag,A,B,C]=odeCody(var_param,var_ind,const_param,const_ind,offset,gain,input,time,exp_data,x0,testflag,w,ntrials)
        %[x,resnorm,residual,exitflag,output]=lsqnonlin(fun,x0,lb,ub,options)
        [bparams1(i,:),bcost1(i),~,exitflag1(i),output1{i}] = ...
            lsqnonlin(@(p)odeChen02(p,var_ind,const_param,const_ind,offset,gain,Input,time,exp_data,[],'none',w,ntrials),...
            startpoints(i,:),plower2,pupper2,options);
    end
    case 'Peng96'
    parfor i=1:montenum
        %[fun,y,flag,A,B,C]=odeCody(var_param,var_ind,const_param,const_ind,offset,gain,input,time,exp_data,x0,testflag,w,ntrials)
        %[x,resnorm,residual,exitflag,output]=lsqnonlin(fun,x0,lb,ub,options)
        [bparams1(i,:),bcost1(i),~,exitflag1(i),output1{i}] = ...
            lsqnonlin(@(p)odePeng96(p,var_ind,const_param,const_ind,offset,gain,Input,time,exp_data,[],'none',w,ntrials),...
            startpoints(i,:),plower2,pupper2,options);
    end
    case 'Cody2'
    parfor i=1:montenum
        %[fun,y,flag,A,B,C]=odeCody(var_param,var_ind,const_param,const_ind,offset,gain,input,time,exp_data,x0,testflag,w,ntrials)
        %[x,resnorm,residual,exitflag,output]=lsqnonlin(fun,x0,lb,ub,options)
        [bparams1(i,:),bcost1(i),~,exitflag1(i),output1{i}] = ...
            lsqnonlin(@(p)odeCody2(p,var_ind,const_param,const_ind,offset,gain,Input,time,exp_data,[],'none',w,ntrials),...
            startpoints(i,:),plower2,pupper2,options);
    end
    end
    %Find the best final cost value out of the montenum starting points.
    %The parameters associated with this cost are the "best" parameters for
    %that trial.
    [bcost,bestind]=min(bcost1);
    bparams=bparams1(bestind,:);
    exitflag=exitflag1(bestind);
    output0=output1{bestind};
  
    switch models{modelind}
        %since the 3 inputs are the same, we simulate y_angles once by
        %setting ntrials=1, not 3. We use the first input vector.
        case 'Cody'
    [~,y_angles]=odeCody(bparams,var_ind,const_param,const_ind,offset,...
        gain,Input(:,1),time,[],[],'simulate',[],1);
        case 'Chen02'
    [~,y_angles]=odeChen02(bparams,var_ind,const_param,const_ind,offset,...
        gain,Input(:,1),time,[],[],'simulate',[],1);
        case 'Peng96'
    [~,y_angles]=odePeng96(bparams,var_ind,const_param,const_ind,offset,...
        gain,Input(:,1),time,[],[],'simulate',[],1);
        case 'Cody2'
    [~,y_angles]=odeCody2(bparams,var_ind,const_param,const_ind,offset,...
        gain,Input(:,1),time,[],[],'simulate',[],1);
    end
    switch models{modelind}
        case {'Cody','Cody2','Chen02'}
    save(strcat(cd2,models{modelind}),'bparams','offset','gain','y_angles','bcost','exitflag','output0','Input','Output','time');
        case 'Peng96'
    save(strcat(strcat(cd2,models{modelind}),num2str(var_ind)),'bparams','offset','gain','y_angles','bcost','exitflag','output0','Input','Output','time','var_ind','const_param','const_ind');        
    end
    cd ..
    end
    cd ..
end
cd ..
end
toc
% delete(gcp('nocreate'))%delete the current parallel pool