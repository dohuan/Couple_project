%This function penalizes y_sim-y_exp across all the ntrials of a given task
%controller is a PID
function [fun,y,flag,A,B,C]=...
    odeCody(var_param,var_ind,const_param,const_ind,offset,gain,input,time,exp_data,x0,testflag,w,ntrials)
%This is a multi-purpose function that takes unknown parameters and known
%parameters and builds the closed-loop head/neck system.  Depending on
%how "testflag" is set, the function either checks for stability of the
%closed-loop A, or simulates the closed-loop system and compares the
%response to experimental data.

%The parameters are
%[delay,Kp,Ki,Kd,k,b,I,tc]
var_param=var_param.*gain+offset;
delay=var_param(var_ind==1);
if isempty(delay), delay=const_param(const_ind==1); if isempty(delay), error('delay not found'); end, end
Kp=var_param(var_ind==2);
if isempty(Kp), Kp=const_param(const_ind==2); if isempty(Kp), error('Kp not found'); end, end
Ki=var_param(var_ind==3);
if isempty(Ki), Ki=const_param(const_ind==3); if isempty(Ki), error('Ki not found'); end, end
Kd=var_param(var_ind==4);
if isempty(Kd), Kd=const_param(const_ind==4); if isempty(Kd), error('Kd not found'); end, end
k=var_param(var_ind==5);
if isempty(k), k=const_param(const_ind==5); if isempty(k), error('k not found'); end, end
b=var_param(var_ind==6);
if isempty(b), b=const_param(const_ind==6); if isempty(b), error('b not found'); end, end
I=var_param(var_ind==7);
if isempty(I), I=const_param(const_ind==7); if isempty(I), error('I not found'); end, end
tc=var_param(var_ind==8);
if isempty(tc), tc=const_param(const_ind==8); if isempty(tc), error('tc not found'); end, end

%Build the plant
Plant=ss([0 1;-k/I -b/I],[0;1/I],[1 0],[0]);
Plant.inputname={'T'};
Plant.outputname={'y'};

%Build the fixed-structure feedback controller
Ksys=tf([Kd Kp Ki],[1 0]);
Ksys.inputname={'e'};
Ksys.outputname='u';

%Input delay
delaysys=tf(1,1,'InputDelay',delay);
% [n,d]=pade(delay,5);
% delaysys=tf(n,d);
delaysys.inputname='r';
delaysys.outputname='rd';

%1st-order muscle dynamics
mdyn=tf(1,[tc 1]);
mdyn.inputname='u';
mdyn.outputname='T';

%Create the summing junction
Sum = sumblk('e = rd - y');

sys0=connect(Plant,Ksys,delaysys,mdyn,Sum,'r',{'y'});
%The three input arguments following sys specify the approximation orders of any input, output, and internal delays of sys, respectively. inf specifies that a delay is not to be approximated.
sys1=pade(sys0,inf,inf,5);%sys1 is a pade-approximated sys0
Ts=time(2)-time(1);
sys1d=c2d(sys1,Ts);
[A,B,C,~]=ssdata(sys1d);

if strcmp(testflag,'buildsys'), fun=[]; y=[]; flag=[]; return; end
%testflag is used to test for stability of the closed-loop system.  If it
%is not set, then we just simulate the system.
if strcmp(testflag,'stability')
    if isstable(sys1)%logical value
        flag='stable';
    else
        flag='unstable';
    end
    y=[];
    fun=[];
    return
end

y=zeros(size(exp_data));
for i=1:ntrials
    [y(:,i),~]=lsim(sys0,input(:,i),time);
end
if strcmp(testflag,'simulate')
    fun=[]; flag=[];
    return;
end
fun=y-exp_data;
fun=fun(:);