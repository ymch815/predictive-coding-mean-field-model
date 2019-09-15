% this function plot the array either by total populations or by
% subpopulations
function pop_plot(xarr,arrs,grpsize,type)
figure; hold on; 
if strcmp(type,'total') || strcmp(type,'subp')
    arr_n = arrs{1}; arr_c = arrs{2};
[arr_np,arr_nerr]=avewithgroup(arr_n,grpsize);
[arr_cp,arr_cerr]=avewithgroup(arr_c,grpsize);
if strcmp(type,'total')
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_np(:,1),arr_nerr(:,1),'linewidth',2)
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_cp(:,1),arr_cerr(:,1),'linewidth',2)
elseif strcmp(type,'subp')
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_np(:,2),arr_nerr(:,2),'linewidth',2)
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_cp(:,2),arr_cerr(:,2),'linewidth',2)
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_np(:,3),arr_nerr(:,3),'linewidth',2)
errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_cp(:,3),arr_cerr(:,3),'linewidth',2)
end
elseif strcmp(type,'other')
    for i=1:size(arrs,2)
        [arr_p,arr_err]=avewithgroup(arrs{i},grpsize);
        errorbar(xarr(ceil(grpsize/2):grpsize:end),arr_p(:,1),arr_err(:,1),'linewidth',2)
    end
        
end
end
