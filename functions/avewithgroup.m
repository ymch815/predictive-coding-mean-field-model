% this function help us to average every N data in an array
function [ave,stda]=avewithgroup(arr,grpsize)
len=size(arr,1);
ave=zeros(len/grpsize,size(arr,2));
stda=zeros(len/grpsize,size(arr,2));

for i=1:len/grpsize
    ave(i,:)=mean(arr(1+(i-1)*grpsize:i*grpsize,:),1);
    stda(i,:)=std(arr(1+(i-1)*grpsize:i*grpsize,:),1)/sqrt(grpsize);
end
end