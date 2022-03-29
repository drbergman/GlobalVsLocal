function [means,stds,means_by_type,cod_intensities] = ...
    codensityFunction(locs,n,types,num_types,cod_edges,cod_intensity_mat,track_intensities)

d = locs-reshape(locs',1,3,[]);
d = squeeze(sqrt(sum(d.^2,2)));
d_sorted = sort(d);
if n<=size(d_sorted,1)
    cod = d_sorted(n,:);
else
    cod = d_sorted(end,:);
end

cod_intensities = zeros(length(cod_edges),num_types+1,num_types);
for ti = num_types:-1:1
    cod_temp = cod(types==ti);
    means(ti) = mean(cod_temp);
    stds(ti) = std(cod_temp);
    if track_intensities
        counts = histcounts(cod_temp,cod_edges,'Normalization','countdensity');
        relfreqs = counts/max(counts);
        cod_intensities(:,4,ti) = cod_intensity_mat*relfreqs';
    end
end

means_by_type = zeros(num_types);
for ti = num_types:-1:1
    d_temp = d(:,types==ti);

    for oi = num_types:-1:1

        d_temp_temp = d_temp(types==oi,:);
        d_temp_temp = sort(d_temp_temp);
        if n<=size(d_temp_temp,1)
            cod = d_temp_temp(n,:);
        else
            cod = d_temp_temp(end,:);
        end
        means_by_type(oi,ti) = mean(cod);

        if track_intensities
            counts = histcounts(cod,cod_edges,'Normalization','countdensity');
            relfreqs = counts/max(counts);
            cod_intensities(:,oi,ti) = cod_intensity_mat*relfreqs';
        end
    end
end