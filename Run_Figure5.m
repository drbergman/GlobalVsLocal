clearvars;

% this will run the simulations that make up Figure 5, exploring the
% agreement between the methods in the IL6 example

%% set up cohort
nsamps = 20;
method = ["local","global"];
base_name = "Figure5_";
min_parfor = 4; % if the timing of these runs is not important, allow for use of a parallel pool if desired
save_output = false; % set to true if you want to save the output

%% set up common parameters
input = basePars_IL6();

%% some parameters you may want to change for this cohort

input.simpars.censor_date = 50 * 24 * 60; % simulate for 50 days (1 day = 24 * 60 minutes)
% input.pars.desired_dt = 1;
% input.pars.TME_size = [20,20,20];
input.initialization_pars.N0 = [100,10,20];

%% run cohort

sz = [length(method),nsamps];
total_runs = prod(sz);
F(1:total_runs) = parallel.FevalFuture;

for mi = length(method):-1:1
    input.method = method(mi);
    for si = nsamps:-1:1
        if total_runs>=min_parfor
            ind = sub2ind(sz,mi,si);
            F(ind) = parfeval(@simPatient_IL6,1,input);
        else
            out(si,mi) = simPatient_IL6(input);
        end
    end
end

if total_runs>=min_parfor
    for i = 1:total_runs
        [idx,next_out] = fetchNext(F);
        [mi,si] = ind2sub(sz,idx);
        out(si,mi) = next_out;
    end
end

%% save output
if save_output
    if ~exist('./data','dir') % make sure the data directory exists
        mkdir('./data')
    end
    filename = nextFileName('data',base_name,3);

    clear F next_out % no need to save these variables

    save(filename,'-v7.3')
end

%% generate panels

nfigs = 1;
figs = gobjects(nfigs,1);
types = ["Stem","Progenitor","TD"];
disp_names = ["Local","Global"];
method_colors = lines(2);
t = out(1,1).tracked.T/(60*24); % time in days
nt = length(t);

%% cell populations
figs(1) = figure("Name","Type Populations");
ax = gobjects(1,3);
y = zeros(nt,nsamps);
ls_type = gobjects(2,3);

for ti = 1:3
    ax(ti) = subplot(1,3,ti);
    hold on
    title(types(ti))
    for mi = 1:length(method)
        for si = 1:nsamps
            y(:,si) = out(si,mi).tracked.NT(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
    end
end
legend(ax(end),ls_type(:,end),'Location','best')
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))

set(ax,'FontSize',20)
for i = 1:numel(ax)
    ax(i).XLabel.String = 'Time (d)';
end

legend(ax(end),ls_type(:,end),'Location','best')
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'FontSize',20)

%% wall time panels
figs(2) = figure("Name","Wall Time Panels");
ax = subplot(1,3,1);
hold on
ylabel('Frequency')
xlabel('Wall Time (h)')
title('Total Time')
t_temp = t;
t_temp(end) = [];
WT = zeros(nsamps,length(method));
ls = gobjects(length(method),1);
for mi = 1:length(method)
    for si = 1:nsamps
        WT(si,mi) = sum(out(si,mi).tracked.wall_time)/3600;
    end
end

for mi = 1:length(method)
    histogram(log10(WT(:,mi)),"FaceColor",method_colors(mi,:),'DisplayName',disp_names(mi))
end
xL = [floor(log10(min(WT,[],'all'))),ceil(log10(max(WT,[],'all')))];
xlim(xL)
xt = ceil(xL(1)):floor(xL(2));
xticks(xt)
xticklabels(10.^xt)

ax(2) = subplot(1,3,2);
hold on
xlabel('Simulation Time (d)')
ylabel(["Cumulative","Wall Time (h)"])
title(["Time on Molecular","Dynamics"])
WT = zeros(nt-1,nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        WT(:,si,mi) = cumsum(out(si,mi).tracked.wall_time_ode+out(si,mi).tracked.wall_time_pde)/3600;
    end
    s = std(WT(:,:,mi),[],2);
    ybar = mean(WT(:,:,mi),2);
    ls(mi) = plot(ax(2),t_temp,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
    patch(ax(2),[t_temp;flip(t_temp)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
end

set(ax(2),'XLim',[0,t_temp(end)])
set(ax(2),'XTick',0:10:input.simpars.censor_date/(60*24))
set(ax(2),'XLim',[0 input.simpars.censor_date/(60*24)])

ax(3) = subplot(1,3,3);
hold on
xlabel("Cells")
ylabel(["Wall Time","Per Update (s)"])
title(["Time by","Number of Agents"])
N = zeros(nt-1,nsamps,length(method));
WT = zeros(nt-1,nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        N(:,si,mi) = sum(out(si,mi).tracked.NT(1:end-1,:),2);
        WT(:,si,mi) = out(si,mi).tracked.wall_time;
    end
    scatter(ax(3),N(:,:,mi),WT(:,si,mi),5,method_colors(mi,:),'Marker','.')
    ls(mi) = plot(ax(3),t_temp,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end
set(ax(3),'XLim',[0,max(N,[],'all')])

%% mean normalized distance to blood vessels
figs(3) = figure("Name","Average Normalized Distance to Blood Vessel");
ax = gobjects(1,3);
y = zeros(nt,nsamps);
for ti = 1:3
    ax(ti) = subplot(1,3,ti);
    hold on
    title(types(ti))
    for mi = 1:length(method)
        for si = 1:nsamps
            y(:,si) = out(si,mi).tracked.mean_z(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
   end
end
legend(ax(end),ls_type(:,end),'Location','best')
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'FontSize',20)
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))

%% il6 concentrations by distance from blood vessel and time
figs(4) = figure("Name","Average IL-6 by Z");
ax = gobjects(3,1);
Z = cell(1,2);
ncolors = 15;
n_percent_colors = 100; % number of colors to use in + %s (same for - %s)
for mi = 1:2
    ax(mi) = subplot(1,3,mi);
    for si = 1:nsamps
        Z{mi}(:,:,si) = out(si,mi).tracked.substrate_by_z(:,:,1)';
    end

    contourf(out(si,mi).tracked.T/(60*24),linspace(0,1,out(si,mi).grid.size(3)),mean(Z{mi},3),"LineColor","none")
    c=colorbar;
    c.Label.String = 'IL-6 Concentration (nM)';
    c.Location = "southoutside";
    title(ax(mi),disp_names(mi))
    if mi==1
        ylabel(["Normalized Distance to","Blood Vessel"])
    end

end

cmap = flipud(bone(ncolors));
colormap(cmap)

ax(end) = subplot(1,3,3);
percent_diffs = 100*(mean(Z{2},3)./mean(Z{1},3)-1);
percent_diffs(isnan(percent_diffs)) = 0;
max_percent_diff = max(abs(percent_diffs),[],'all');
[bounds_percent_diff(1),bounds_percent_diff(2)] = bounds(percent_diffs(:));

if isequal(bounds_percent_diff,[0,0])
    bounds_percent_diff = [-1 1];
end
cmap_pd_pos = bone(n_percent_colors);
cmap_pd_neg = bone(n_percent_colors);

if abs(bounds_percent_diff(1))>abs(bounds_percent_diff(2)) % then the lower bound is negative and larger than the upper bound
    upper_ind = floor((1-abs(bounds_percent_diff(2))/abs(bounds_percent_diff(1)))*n_percent_colors);
    cmap_pd_pos(1:upper_ind,:) = [];
else
    upper_ind = floor((1-abs(bounds_percent_diff(1))/abs(bounds_percent_diff(2)))*n_percent_colors);
    cmap_pd_neg(1:upper_ind,:) = [];
end

cmap_pd = cat(1,cmap_pd_neg,...
    flipud(cmap_pd_pos));
contourf(out(si,mi).tracked.T/(60*24),linspace(0,1,out(si,mi).grid.size(3)),percent_diffs,n_percent_colors,"LineColor","none")
ax(end).Children.LevelList = linspace(-100,100,2*n_percent_colors+1);
colormap(ax(end),cmap_pd)
caxis(ax(end),bounds_percent_diff)
c=colorbar;
c.Label.String = 'Percent Difference';
c.Location = "southoutside";
title(["Global Relative","to Local"])

for i = 1:numel(ax)
    ax(i).XLabel.String = 'Time (d)';
end
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))


