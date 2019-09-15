% this function define what we do in a trial
% it describes the key features of our model

function [xarr,returnpara,dynamics] = trials(state, trial, mode, returntype, ...
    plotma, plotma_1, itr, itrmin, itrmax, rep, noise)

global base T dt pre N

% Set stimulus period
if strcmp('evoked',trial) || strcmp('spontaneous',trial)
t1 = 400+base; 
t2 = 600+base; 
% t1 and t2 are the start and ending time of stimulus
elseif strcmp('placebo',trial)
t1 = 200+base; 
t2 = 400+base;
elseif strcmp('prediction',trial(1:10))
t1 = 400+base; 
t2 = 600+base;
end

% Model parameters
wEE = 22; 
wEI = 22;

if strcmp('chronic',state)
long_range = 0.2; % >0.1
per = 0.3;  % 0.3, the percentage of neurons in ACC which directly receive the input from S1
prS1= 0.6; % 0.52
prNS1= 0.2; % 0.20
pr= 0.25;
end

if strcmp('naive',state)
long_range = 0.1; % 0.1
per = 0.2;  % 0.2
prS1= 0.35; % 0.35
prNS1=0.14; % 0.14
pr=0.1;
end

gS1=2;
gACC=5;

mysize = [2,-1.5,per,(1-per)]; % S1/ACC(alpha), I/E(rho), pain-responsive-S1/all-pain-responsive(per)

% set conneciton matrix
w = zeros(5,5);
w(1,1)=wEE; w(2,2)=wEE*mysize(2); w(3,3)=wEE*mysize(3)/mysize(1); 
w(4,4)=wEE*mysize(4)/mysize(1); w(5,5)=wEE*mysize(2)/mysize(1);
w(1,2)=wEI; w(2,1)=wEI*mysize(2);
w(1,3)=wEE*long_range; w(3,1)=wEE*long_range;
if strcmp('no feedback',mode)
w(3,1)=0;
end
w(3,4)=wEE*mysize(3)/mysize(1); w(4,3)=wEE*mysize(4)/mysize(1);
w(3,5)=wEI*mysize(3)/mysize(1); w(5,3)=wEI*mysize(2)/mysize(1);
w(4,5)=wEI*mysize(4)/mysize(1); w(5,4)=wEI*mysize(2)/mysize(1);

sigmaE = [0.5,0.7]; % 0.6 0.5naive
sigmaI = [0.5,0.7];
hE = [4,3]; % 4naive
hI = [4,3];

epsilonE = 0.005; 
epsilonI = 0.005; 
if noise==0
epsilonE = 0.00; 
epsilonI = 0.00; 
end
epsilonZ = 0.1; % 0.1

tauSE = 3;
tauSI = 10;
satur = 4; % 5, 3.5, 6naive

tauE1 = 1; %0.2 4
tauI1 = 3;  % 1.5 12
tauE21 = 3; %7 17/3
tauE22 = 3; %7 17/3
tauI2 = 18; % 14 7

delayt = 20; % delay from S1 to ACC in unit of ms
delay = delayt/dt;

delayxt = 20; % delay from x to z
delayx = delayxt/dt;

delayunit=30; % delay unit in ms
delayfb = zeros(N+1,1);
delaygenerator = lognrnd(log(80),0.6,ceil((N+1)/(delayunit/dt)),1);
for i=0:size(delaygenerator,1)-1
delayfb(i*delayunit/dt+1:min((i+1)*delayunit/dt,N+1))=delaygenerator(i+1);
end

thresd=0.91;
itrarr=linspace(itrmin,itrmax,itr);

xarr_ = zeros(itr*rep,1);
Recell = cell([1 size(returntype,2)]);
for i=1:size(returntype,2)
    Recell{i}=zeros(itr*rep,3);
end


for itrt=0:size(itrarr,2)-1
    setvalue = itrarr(itrt+1); 
for repeat=1:rep
        
% Set vectors
E1 = zeros(N+1,1);
I1 = zeros(N+1,1);
E21 = zeros(N+1,1);
E22 = zeros(N+1,1);
I2 = zeros(N+1,1);
t = linspace(0,T,N+1);
S = zeros(N+1,5);
x = zeros(N+1,1);

if strcmp('evoked', trial)
z = zeros(N+1,1);
tpoint=t1/dt; %the staring point of z integration
zstart=0; % the starting point of z undating
t7=t1; % starting point of peak integration
x(t1/dt:t2/dt) = setvalue; % Set stimulus
paraA = 2000; 
Zth=300; 
elseif strcmp('spontaneous',trial)
z = zeros(N+1,1);
z(base/dt+1,1)=setvalue;
zstart=base;
tpoint=base/dt+1; 
t7=base; % starting point of peak integration
paraA = 2000; 
Zth=400;
elseif strcmp('placebo',trial)
z = zeros(N+1,1);
z(base/dt+1,1)=setvalue;
zstart=base;
tpoint=t1/dt;
t7=t1; % starting point of peak integration
x(t1/dt:t2/dt) = 2.5;
paraA = 2000;
Zth= 500; 
elseif strcmp('prediction',trial(1:10))
z = zeros(N+1,1);
zstart=t1;
tpoint=t1/dt;
t7=t1; % starting point of peak integration
paraA = 2000; 
Zth=400; 
if strcmp('prediction_precise',trial)
z(base/dt+1:t1/dt+1,1)=setvalue;
x(t1/dt:t2/dt) = setvalue; 
elseif strcmp('prediction_wrong',trial)
z(base/dt+1:t1/dt+1,1)=0;
x(t1/dt:t2/dt) = setvalue;
elseif strcmp('prediction_const',trial(1:16))
z(base/dt+1:t1/dt+1,1)=setvalue;
x(t1/dt:t2/dt) = str2double(trial(18:end)); 
end
end

p_reset=0; % initialize the point of reset

% Initialize Gaussian white noise
if epsilonE == 0 && epsilonI == 0 
  noise = zeros(5,N);
else
  noise = normrnd(0,1,[5,N]);
end

f = @(x,sigma,h) 1./(1+exp(-sigma*(x-h)));

% compute latent variable z and the point of reset
tauz = paraA./(1+exp(x));
noisez = normrnd(0,1,[1,N])*epsilonZ;
count = 0;
for i=zstart/dt+1:N
   if i>delayx
    z(i+1) = z(i) + dt*(-z(i)+x(i-delayx)+noisez(i)/sqrt(dt))/tauz(i);
    else
    z(i+1) = z(i) + dt*(-z(i)+noisez(i)/sqrt(dt))/tauz(i);
    end
    if i+1>tpoint && sum(z(tpoint:i+1))*dt>Zth
        z(i+1) = 0;
        tpoint = i+1;
        count=count+1;
        if z(i)>0.3 || count==1 % ignore subsequent rise that is too small
        p_reset = i+1; 
        x(p_reset-delayx:end)=0; % set x to 0 if withdraw
        end
    end
end

% put x and z as inputs to the network
PE1 = abs(x-z)*gS1; 
PE21 = z*gACC*prS1;
PE22 = z*gACC*prNS1;
PI1 = PE1;
PI2 = z*gACC*pr;

for i=1:N
    % Iterate using Euler method
    % Equations for population 1
    if strcmp('feedback',mode)
        if i<=delayfb(i)/dt
        FE1 = f(w(:,1)'*S(i,:)'-w(3,1)*S(i,3)+PE1(i),sigmaE(1),hE(1));
        else
        FE1 = f(w(:,1)'*S(i,:)'-w(3,1)*S(i,3)+w(3,1)*S(ceil(i-delayfb(i)/dt),3)+PE1(i),sigmaE(1),hE(1));
        end
    else
        FE1 = f(w(:,1)'*S(i,:)'+PE1(i),sigmaE(1),hE(1));
    end
    
    FI1 = f(w(:,2)'*S(i,:)'+PI1(i),sigmaI(1),hI(1));
  
    S(i+1,1) = S(i,1) + dt*(-S(i,1)+satur*(1-S(i,1))*E1(i)...
        +epsilonE/sqrt(dt)*noise(1,i))/tauSE;
    S(i+1,2) = S(i,2) + dt*(-S(i,2)+satur*(1-S(i,2))*I1(i)...
        +epsilonI/sqrt(dt)*noise(2,i))/tauSI;
   
    E1(i+1) = E1(i) + dt*(-E1(i)+FE1)/tauE1;
    I1(i+1) = I1(i) + dt*(-I1(i)+FI1)/tauI1;
   
    % Equations for population 2
    if i<=delay
    FE21 = f(w(:,3)'*S(i,:)'-w(1,3)*S(i,1)+PE21(i),sigmaE(2),hE(2));
    else
    FE21 = f(w(:,3)'*S(i,:)'-w(1,3)*S(i,1)+PE21(i)+...
        w(1,3)*S(i-delay,1),sigmaE(2),hE(2));
    end
    
    FE22 = f(w(:,4)'*S(i,:)'+PE22(i),sigmaE(2),hE(2));

    FI2 = f(w(:,5)'*S(i,:)'+PI2(i),sigmaI(2),hI(2));
   
    S(i+1,3) = S(i,3) + dt*(-S(i,3)+satur*(1-S(i,3))*E21(i)...
        +epsilonE/sqrt(dt)*noise(3,i))/tauSE;
    S(i+1,4) = S(i,4) + dt*(-S(i,4)+satur*(1-S(i,4))*E22(i)...
        +epsilonI/sqrt(dt)*noise(4,i))/tauSE;
    S(i+1,5) = S(i,5) + dt*(-S(i,5)+satur*(1-S(i,5))*I2(i)...
        +epsilonI/sqrt(dt)*noise(5,i))/tauSI;
   
    E21(i+1) = E21(i) + dt*(-E21(i)+FE21)/(tauE21);
    E22(i+1) = E22(i) + dt*(-E22(i)+FE22)/(tauE22);
    I2(i+1) = I2(i) + dt*(-I2(i)+FI2)/(tauI2);
    
end

E21r = E21; 
E1 = S(:,1);
E21 = S(:,3);
E22 = S(:,4);

method='pchip';
E1e = myenvelope(E1,t/dt,[0.2,0.8],100,ceil(size(E1,1)/350),method);
E21e = myenvelope(E21,t/dt,[0.3,0.7],300,ceil(size(E21,1)/100),method);
E21re = myenvelope(E21r,t/dt,[0.3,0.7],300,ceil(size(E21r,1)/100),method);
E22e = myenvelope(E22,t/dt,[0.3,0.7],300,ceil(size(E22,1)/100),method);

E1m = (E1e(1,:)+E1e(2,:))/2;
E21m = (E21e(1,:)+E21e(2,:))/2;
E21rm = (E21re(1,:)+E21re(2,:))/2;
E22m = (E22e(1,:)+E22e(2,:))/2;

t3=1+base; % the starting point of preS1 integration
t4=p_reset*dt; % the end point of preS1 integration
timelen=t4/dt-t3/dt+1; % #points of preS1 integration
timelen_acc=1000/dt; % #points of postACC integration
t5=500; % starting point of baseline calculation
t6=1500; % ending point of baseline calculation
timelen_base=t6/dt-t5+1; % # points of baseline calculation
t8=p_reset*dt; % ending point of peak integration
peaklen=t8/dt-t7/dt+1; % # points of peak integration
t9=1+base; % the starting point of beforeS1 integration
t10=t1; % the end point of beforeS1 integration
timelen_bef=t10/dt-t9/dt+1; % #points of beforeS1 integration

peakposit21=find(E21(p_reset-40/dt:p_reset)==max(E21(p_reset-40/dt:p_reset)));
peakposit22=find(E22(p_reset-40/dt:p_reset)==max(E22(p_reset-40/dt:p_reset)));
peakposit21=p_reset-40/dt+peakposit21-1;
peakposit22=p_reset-40/dt+peakposit22-1;

% update values
xarr_(itrt*rep+repeat) = setvalue; 

for j=1:size(returntype,2)
    type = returntype{j};
switch type
    case 'preS1'
        Recell{j}(itrt*rep+repeat,1)=sum(E1m(t3/dt:t4/dt))/timelen;
    case 'beforeS1'
        Recell{j}(itrt*rep+repeat,1)=sum(E1m(t9/dt:t10/dt))/timelen_bef;
    case 'duringS1'
        Recell{j}(itrt*rep+repeat,1)=sum(E1m(t7/dt:t8/dt))/peaklen;
    case 'postACC'
        post_acc_dir = sum(E21m(peakposit21:peakposit21+timelen_acc-1))/timelen_acc; 
        post_acc_ind = sum(E22m(peakposit22:peakposit22+timelen_acc-1))/timelen_acc;
        post_acc = post_acc_dir+post_acc_ind;
        Recell{j}(itrt*rep+repeat,1)=post_acc;
        Recell{j}(itrt*rep+repeat,2)=post_acc_dir;
        Recell{j}(itrt*rep+repeat,3)=post_acc_ind;
    case 'peakACC'
        peak_acc_dir = sum(E21m(t7/dt:t8/dt))/peaklen; 
        peak_acc_ind = sum(E22m(t7/dt:t8/dt))/peaklen; 
        peak_acc = peak_acc_dir+peak_acc_ind;
        Recell{j}(itrt*rep+repeat,1)=peak_acc;
        Recell{j}(itrt*rep+repeat,2)=peak_acc_dir;
        Recell{j}(itrt*rep+repeat,3)=peak_acc_ind;
    case 'sustain'
        Recell{j}(itrt*rep+repeat,1)=size(find((E21m(1+base/dt:p_reset)+...
            E22m(1+base/dt:p_reset))>thresd),2)/...
           size((E21m(1+base/dt:p_reset)+E22m(1+base/dt:p_reset)),2);
    case 'maximum'
        [Recell{j}(itrt*rep+repeat,2),~] = max(E21m(t7/dt:t8/dt));
        [Recell{j}(itrt*rep+repeat,3),~] = max(E22m(t7/dt:t8/dt));
        [Recell{j}(itrt*rep+repeat,1),~] = max(E21m(t7/dt:t8/dt)+E22m(t7/dt:t8/dt));
    case 'latency'
        [~,latency_acc_dir] = max(E21(t7/dt:t8/dt));
        [~,latency_acc_ind] = max(E22(t7/dt:t8/dt));
        [~,latency_acc] = max(E21(t7/dt:t8/dt)+E22(t7/dt:t8/dt));
        Recell{j}(itrt*rep+repeat,1)=latency_acc*dt;
        Recell{j}(itrt*rep+repeat,2)=latency_acc_dir*dt;
        Recell{j}(itrt*rep+repeat,3)=latency_acc_ind*dt;
    case 'baselineACC'
        acc_S1base = sum(E21m(t5/dt:t6/dt))/timelen_base;
        acc_nS1base = sum(E22m(t5/dt:t6/dt))/timelen_base;
        acc_base = acc_S1base+acc_nS1base;
        Recell{j}(itrt*rep+repeat,1)=acc_base;
        Recell{j}(itrt*rep+repeat,2)=acc_S1base;
        Recell{j}(itrt*rep+repeat,3)=acc_nS1base;
    otherwise
end
end

if itrt==0
    fprintf('"%s" condition under "%s" pain, "%s"\n',trial, state, mode);
end
fprintf('x(or z) magnitude=%4.2f,repeat %2.0f time(s)\n',...
    setvalue,repeat);
if itrt==size(itrarr,2)-1 && itr*rep>1
    fprintf('"%s" condition under "%s" pain, "%s"\n',trial, state, mode);
end

if plotma_1==1
figure
subplot(2,1,1)
hold on
plot(t(pre),E21r(pre),'r-','LineWidth',1);
plot(t(pre),E21(pre),'b-','LineWidth',1);
plot(t(pre),E21re(:,pre),'r--','LineWidth',1);
plot(t(pre),E21rm(pre),'r-','LineWidth',2);
plot(t(pre),E21e(:,pre),'b--','LineWidth',1);
plot(t(pre),E21m(pre),'b-','LineWidth',2);
myplot('Time (ms)', {'Original signal','(a.u.)'}, 20,...
    {'Firing rate','Synaptic activ.'}, 20,'Northeast',...
    [pre(1)*dt pre(end)*dt 0 1.2])
set(legend,'NumColumns',2);set(gca,'xticklabel',{[]});

subplot(2,1,2)
hold on
plot(t(pre),E21rm(pre),'r-','LineWidth',2);
plot(t(pre),E21m(pre),'b-','LineWidth',2);
myplot('Time (ms)', {'Middle line','(a.u.)'}, 20,...
    {'Firing rate','Synaptic activ.'}, 20,'Northeast',...
    [pre(1)*dt pre(end)*dt 0 0.8])
set(legend,'NumColumns',2);

figure
hold on 
scatter(E21r(pre(1):30:pre(end)),E21(pre(1):30:pre(end)),60,'filled','MarkerFaceColor',[0.8 0 0.5])
box on; grid off; xlabel('Firing rate variable (a.u.)');
ylabel({'Synaptic activ. variable (a.u.)'});
set(gca,'fontsize',20);set(gca,'linewidth',2);
corr_coe=corr(E21r(pre(1):30:pre(end)), E21(pre(1):30:pre(end)),'type','kendall');
fprintf("the correlation is: %4.2f \n",corr_coe)
end

if plotma==1
figure
subplot(5,1,1)
hold on
plot(t(pre),x(pre),'LineWidth',2);
plot(t(pre),z(pre),'LineWidth',2);
myplot('', 'Amp (a.u.)', 12, {'x','z'}, 10,'Northeast',[pre(1)*dt pre(end)*dt ...
    min(min(z(pre)),min(x(pre)))-0.2 max(max(z(pre)),max(x(pre)))]+0.2)
set(gca,'xticklabel',{[]});

subplot(5,1,2)
hold on
plot(t(pre),E1(pre),'k-','LineWidth',1);
plot(t(pre),E1e(:,pre),'k-','LineWidth',1);
plot(t(pre),E1m(:,pre),'k-','LineWidth',2);
myplot('', '', 12,{'S1-E'}, 10,'Northeast',[pre(1)*dt pre(end)*dt 0 1])
set(gca,'xticklabel',{[]});

subplot(5,1,3)
hold on
plot(t(pre),E21(pre),'r-','LineWidth',1);
plot(t(pre),E21e(:,pre),'r-','LineWidth',1);
plot(t(pre),E21m(:,pre),'r-','LineWidth',2);
myplot('', 'Synaptic activity s (a.u.)', 12,{'ACC-E (S1 input)'}, 10,'Northeast',[pre(1)*dt pre(end)*dt 0 1])
set(gca,'xticklabel',{[]});

subplot(5,1,4)
hold on
plot(t(pre),E22(pre),'b-','LineWidth',1);
plot(t(pre),E22e(:,pre),'b-','LineWidth',1);
plot(t(pre),E22m(:,pre),'b-','LineWidth',2);
myplot('', '', 12,{'ACC-E (no S1 input)'}, 10,'Northeast',[pre(1)*dt pre(end)*dt 0 1])
set(gca,'xticklabel',{[]});

subplot(5,1,5)
hold on
plot(t(pre),E1m(:,pre),'k-','LineWidth',2);
plot(t(pre),E21m(:,pre),'r-','LineWidth',2);
plot(t(pre),E22m(:,pre),'b-','LineWidth',2);
scatter(peakposit21*dt,0,50,'k*');
scatter(peakposit21*dt+timelen_acc*dt,0,50,'k*');
myplot('Time(ms)', '', 12, {'','Middle Lines',''}, 10,'Northeast',...
    [pre(1)*dt pre(end)*dt 0 max([max(E1m(:,pre)) max(E21m(:,pre)) max(E22m(:,pre))])+0.1])
set(legend,'NumColumns',1);
end

end
end
returnpara=Recell;
xarr=xarr_;
dynamics = E21m+E22m; 
end
