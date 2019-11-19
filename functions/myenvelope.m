function Evl = myenvelope(Sig,x,maxh,bet,res,method)
% Sig: signal array; x: x-axis array; maxh: min peak value
% bet: min distance between 2 peaks; res: timepoints for each step
% method: method of interpolate
% Evl: high and low envelope of the signal
len=size(Sig,1);
itr = ceil(len/res);
high=[]; idxhigh=[];
low=[]; idxlow=[];
for i=0:itr-1
range=i*res+1:min(i*res+res,len);
highh=maxh(2)*max(Sig(range))+(1-maxh(2))*min(Sig(range));
lowh=maxh(1)*max(Sig(range))+(1-maxh(1))*min(Sig(range));
[hight,idxhight] = findpeaks(Sig(range),'MinPeakHeight',highh);
[lowt,idxlowt] = findpeaks(-Sig(range),'MinPeakHeight',-lowh);
low=[low; -lowt]; high=[high; hight];
idxlow=[idxlow; idxlowt+i*res]; idxhigh=[idxhigh; idxhight+i*res];
% TF1 = islocalmax(Sig(range), 'FlatSelection', 'first');
% plat=Sig(TF1);
% indplat=find(TF1==1);
% low=[low; plat]; high=[high; plat];
% idxlow=[idxlow; indplat+i*res]; idxhigh=[idxhigh; indplat+i*res];
end

for ii=1:size(idxhigh,1)-1
    if (idxhigh(ii+1)-idxhigh(ii))<bet
        if high(ii+1)>high(ii)
            idxhigh(ii)=0;
        else idxhigh(ii+1)=0;
        end
    end
end
id1 = idxhigh(:,1)==0;
idxhigh(id1,:)=[]; high(id1,:)=[];

for ii=1:size(idxlow,1)-1
    if (idxlow(ii+1)-idxlow(ii))<bet
        if low(ii+1)<low(ii)
            idxlow(ii)=0;
        else idxlow(ii+1)=0;
        end
    end
end
id2 = idxlow(:,1)==0;
idxlow(id2,:)=[]; low(id2,:)=[];

Evlh = interp1(idxhigh,high,x,method);
Evll = interp1(idxlow,low,x,method);
Evl=[Evlh;Evll];