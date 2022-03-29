function [means,stds,means_by_type] = ...
    codensityFunction(locs,n,types,num_types)

d = locs-reshape(locs',1,3,[]);
d = squeeze(sqrt(sum(d.^2,2)));
d_sorted = sort(d);
if n<=size(d_sorted,1)
    cod = d_sorted(n,:);
else
    cod = d_sorted(end,:);
end

for ti = num_types:-1:1
    cod_temp = cod(types==ti);
    means(ti) = mean(cod_temp);
    stds(ti) = std(cod_temp);
end

means_by_type = zeros(num_types);
for ti = num_types:-1:1
    d_temp = d(:,types==ti);
    if isempty(d_temp) % if no cells of this type are present
        means_by_type(:,ti) = NaN;
        continue;
    end

    for oi = num_types:-1:1

        d_temp_temp = d_temp(types==oi,:);
        if isempty(d_temp_temp) % if no cells of this type are present
            means_by_type(oi,ti) = NaN;
            continue;
        end
        d_temp_temp = sort(d_temp_temp);
        if n<=size(d_temp_temp,1)
            cod = d_temp_temp(n,:);
        else
            cod = d_temp_temp(end,:);
        end
        means_by_type(oi,ti) = mean(cod);
    end
end