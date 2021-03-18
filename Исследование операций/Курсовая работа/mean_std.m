function [mean_vv, std_vv] = mean_std(i, vv, depth)
mean_vv = mean(vv(max([1,i - depth]):i));
std_vv = sqrt(var(vv(max([1,i - depth]):i)));
end