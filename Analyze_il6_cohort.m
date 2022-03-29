clearvars;
close all

data_file = "il6_cohort001";
load(sprintf('data/%s.mat',data_file))

types = ["Stem","Progenitor","TD"];
if length(method)==3
    disp_names = ["Local","Global + PDE","Global"];
elseif length(method)==2
    if use_pde_in_global
        disp_names = ["Local","Global + PDE"];
    else
        disp_names = ["Local","Global"];
    end
elseif length(method)==1
    if method=="global"
        disp_names = "Global";
    elseif method=="local"
        disp_names = "Local";
    else
        error("not sure what the methods are")
    end
else
    error("not sure what the methods are")
end

if ~exist("nsamps","var")
    nsamps = size(out,1);
end

method_colors = cubehelix(2,2,-1.5,1,1,[0,.6],[.4,.6]);
type_colors = cubehelix(3,2,-1.5,1,1,[0,.6],[.4,.6]);

events = 0:4;
event_colors = [65,145,79;... % new tumor cells
                121,121,121;... % proliferating cells
                219,59,53;... % apoptotic cells
                95,48,140;... % movement
                0,0,0;... % resting
                68,55,31]/255; 
type_markers = {'o','s','^'};

%%
name = "scatter_plots";
for mi = 1:length(method)
    figure("Name","Scatter Plot")
    ax = gca;
    hold on
    tumors = out(1,mi).tumors;
    I = false(size(tumors,1),3);
    tum_locs = tumors(:,input.inds.location_inds);
    num_prolif_attempts = sum(tumors(:,input.inds.event_ind)==1);
    num_prolifs = sum(tumors(:,input.inds.event_ind)==0);
    [event_color_ind,~,~] = find((tumors(:,input.inds.event_ind)==events)');
    prolif_ind = find(event_color_ind==2);
    prolif_ind(tumors(prolif_ind,input.inds.type_ind)==3) = []; % any TD cells that proliferated must have been progenitor cells and so were not contact inhibited
    contact_inhibited_ind = randsample(prolif_ind,num_prolif_attempts-num_prolifs,false);
    event_color_ind(contact_inhibited_ind) = 6;
    for ti = 1:3
        I(:,ti) = tumors(:,input.inds.type_ind)==ti;
%         scatter3(tum_locs(I(:,ti),1),tum_locs(I(:,ti),2),tum_locs(I(:,ti),3),80,type_colors(ti,:),'filled','Marker',type_markers{ti})
        scatter3(tum_locs(I(:,ti),1),tum_locs(I(:,ti),2),tum_locs(I(:,ti),3),160,event_colors(event_color_ind(I(:,ti)),:),'filled','Marker',type_markers{ti})
    end
    ax.Box = 'off';
    ax.NextPlot = 'add';
    ax.Color = ax.Parent.Color;
    ax.XTick = [];
    ax.YTick = [];
    ax.ZTick = [];
    ax.XColor = 'none';
    ax.YColor = 'none';
    ax.ZColor = 'none';
    view(ax,-34,28)


%     savefig(sprintf('figs/%s%s_%s',name,method(mi),data_file))
%     print(sprintf('figs/jpgs/%s%s_%s',name,method(mi),data_file),'-djpeg')
%     print(sprintf('figs/epss/%s%s_%s',name,method(mi),data_file),'-depsc')

end


%%
name = "type_populations";
figure("Name","Type Populations");
ax = gobjects(length(method),3);
for ti = 1:3
    for mi = 1:length(method)
        ax(mi,ti) = subplot(length(method),3,r2c(length(method),3,[mi,ti]));
        hold on
        if mi==1
            title(types(ti))
        end
        if ti==1
            ylabel(out(1,mi).method,'FontWeight','bold')
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y = out(si,mi).tracked.NT(:,ti);
            plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
        end
        plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
    end
    normalizeYLims(ax)
end
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))
set(ax,'FontSize',20)

% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
name = "type_population_mean";
figure("Name","Type Populations");
ax_type = gobjects(1,3);
ax_method = gobjects(1,length(method));
for ti = 1:3
    ax_type(ti) = subplot(2,3,ti);
    hold on
    title(types(ti))
    for mi = 1:length(method)
        if ti==1
            method_subplot_width = 3/length(method)-1;
            ax_method(mi) = subplot(2,3,4+(mi-1)*(1+method_subplot_width) + [0,method_subplot_width]);
            title(method(mi))
            hold on
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y(:,si) = out(si,mi).tracked.NT(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax_type(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax_type(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
        ls_method(mi,ti) = plot(ax_method(mi),t,ybar,'Color',type_colors(ti,:),'LineWidth',2,'DisplayName',types(ti));
        patch(ax_method(mi),[t;flip(t)],[ybar-s;flip(ybar+s)],type_colors(ti,:),'EdgeColor','none','FaceAlpha',0.6)
    end
    normalizeYLims([ax_type,ax_method])
end
legend(ax_type(end),ls_type(:,end),'Location','best')
set(ax_type,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax_type,'XTick',0:10:input.simpars.censor_date/(60*24))
set(ax_method,'XTick',0:10:input.simpars.censor_date/(60*24))

set(ax_type,'FontSize',20)
for i = 1:numel(ax_type)
    ax_type(i).XLabel.String = 'Time (d)';
end

legend(ax_type(end),ls_type(:,end),'Location','best')
legend(ax_method(end),ls_method(end,:),'Location','best')
set([ax_type,ax_method],'XLim',[0 input.simpars.censor_date/(60*24)])
set([ax_type,ax_method],'FontSize',20)
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

fnew = figure;
fnew.Units = "inches";
res_factor = 2;
fnew.Position = [0 0 4 1.5]*res_factor;
for i = 1:3
    ax(i) = subplot(1,3,i);
    copyobj(ax_type(i).Children,ax(i))
    ax(i).XLabel.String = 'Time (d)';
    title(types(i))
    if i==1
        ylabel('Cells')
    end
end
legend(ax(1),flip(ax(1).Children(2:2:end)),'Location','Best')
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))

normalizeYLims(fnew)
set(ax,'FontSize',6*res_factor)

% savefig(sprintf('figs/%s_%s_justtype',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_justtype',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_justtype',name,data_file),'-depsc')

fnew.Position(4) = 0.5*fnew.Position(4);
% savefig(sprintf('figs/%s_%s_justtype_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_justtype_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_justtype_short',name,data_file),'-depsc')

%%
name = "wall_time";
f = figure("Name","Wall Time");
f.Units = "inches";
res_factor = 2;
f.Position = [0 0 4 1.5]*res_factor;
ax = subplot(1,3,1);
hold on
ylabel('Frequency')
xlabel('Wall Time (h)')
title('Total Time')
WT = zeros(nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        WT(si,mi) = sum(out(si,mi).tracked.wall_time)/3600;
    end
    histogram(WT(:,mi),"FaceColor",method_colors(mi,:),'DisplayName',disp_names(mi))
end
% minx = floor(log10(min(WT,[],'all')));
minx = 0;
max_power_10 = floor(log10(max(WT,[],'all')));
max_power_10 = round(max(WT,[],'all'),-max_power_10);
set(ax,'FontSize',6*res_factor,'XLim',[minx,max_power_10])
legend('Location','NorthEast')

ax = subplot(1,3,2);
hold on
xlabel('ODE Wall Time (h)')
title('Time on ODEs')
WT = zeros(nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        WT(si,mi) = sum(out(si,mi).tracked.wall_time_ode)/3600;
    end
    histogram(WT(:,mi),"FaceColor",method_colors(mi,:),'DisplayName',disp_names(mi))
end
% minx = floor(log10(min(WT,[],'all')));
minx = 0;
max_power_10 = floor(log10(max(WT,[],'all')));
max_power_10 = round(max(WT,[],'all'),-max_power_10);
set(ax,'FontSize',6*res_factor,'XLim',[minx,max_power_10])
% legend('Location','NorthEast')
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

% f = figure("Name","PDE Wall Time");
% f.Units = "inches";
% res_factor = 2;
% f.Position = [0 0 4/3 1.5]*res_factor;
ax = subplot(1,3,3);
hold on
xlabel('PDE Wall Time (h)')
title('Time on PDEs')
WT = zeros(nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        WT(si,mi) = sum(out(si,mi).tracked.wall_time_pde)/3600;
    end
    histogram(WT(:,mi),"FaceColor",method_colors(mi,:),'DisplayName',disp_names(mi))
end
% minx = floor(log10(min(WT,[],'all')));
minx = 0;
max_power_10 = max(-16,floor(log10(max(WT,[],'all'))));
maxx = ceil(max(WT,[],'all')/10^max_power_10);
if maxx==0
    maxx=1;
end
set(ax,'FontSize',6*res_factor,'XLim',[minx,maxx])
% legend('Location','NorthEast')
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

f.Position(4) = 0.5*f.Position(4);
% savefig(sprintf('figs/%s_%s_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_short',name,data_file),'-depsc')

%%
name = "wall_time_series";
f = figure("Name","Wall Time");
f.Units = "inches";
res_factor = 2;
f.Position = [0 0 4 1.5]*res_factor;
ax = subplot(1,3,1);
hold on
ylabel(["Cumulative","Wall Time (h)"])
xlabel('Simulation Time (d)')
title('Total Time')
t = out(1,1).tracked.T/(60*24);
t(end) = [];
WT = zeros(length(t),nsamps,length(method));
ls = gobjects(length(method),1);
for mi = 1:length(method)
    for si = 1:nsamps
        WT(:,si,mi) = cumsum(out(si,mi).tracked.wall_time)/3600;
    end
    s = std(WT(:,:,mi),[],2);
    ybar = mean(WT(:,:,mi),2);
    ls(mi) = plot(ax,t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
    patch(ax,[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
end
set(ax,'FontSize',6*res_factor,'XLim',[0,t(end)])
legend(ls,'Location','NorthEast')

ax = subplot(1,3,2);
hold on
xlabel('Simulation Time (d)')
title('Time on ODEs')
for mi = 1:length(method)
    for si = 1:nsamps
        WT(:,si,mi) = cumsum(out(si,mi).tracked.wall_time_ode)/3600;
    end
    s = std(WT(:,:,mi),[],2);
    ybar = mean(WT(:,:,mi),2);
    ls(mi) = plot(ax,t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
    patch(ax,[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
end

set(ax,'FontSize',6*res_factor,'XLim',[0,t(end)])

ax = subplot(1,3,3);
hold on
xlabel('Simulation Time (d)')
title('Time on PDEs')
for mi = 1:length(method)
    for si = 1:nsamps
        WT(:,si,mi) = cumsum(out(si,mi).tracked.wall_time_pde)/3600;
    end
    s = std(WT(:,:,mi),[],2);
    ybar = mean(WT(:,:,mi),2);
    ls(mi) = plot(ax,t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
    patch(ax,[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
end
set(ax,'FontSize',6*res_factor,'XLim',[0,t(end)])

% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

f.Position(4) = 0.5*f.Position(4);
% savefig(sprintf('figs/%s_%s_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_short',name,data_file),'-depsc')

%%
name = "wall_time_combo";
f = figure("Name","Wall Time");
f.Units = "inches";
res_factor = 2;
f.Position = [0 0 4 1.5]*res_factor;
ax = subplot(1,3,1);
hold on
ylabel('Frequency')
xlabel('Wall Time (h)')
title('Total Time')
t = out(1,1).tracked.T/(60*24);
t(end) = [];
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
set(ax,'FontSize',6*res_factor)

ax = subplot(1,3,2);
hold on
xlabel('Simulation Time (d)')
ylabel(["Cumulative","Wall Time (h)"])
title(["Time on Molecular","Dynamics"])
t = out(1,1).tracked.T/(60*24);
t(end) = [];
WT = zeros(length(t),nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        WT(:,si,mi) = cumsum(out(si,mi).tracked.wall_time_ode+out(si,mi).tracked.wall_time_pde)/3600;
    end
    s = std(WT(:,:,mi),[],2);
    ybar = mean(WT(:,:,mi),2);
    ls(mi) = plot(ax,t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
    patch(ax,[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
end

set(ax,'FontSize',6*res_factor,'XLim',[0,t(end)])
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])

ax = subplot(1,3,3);
hold on
xlabel("Cells")
ylabel(["Wall Time","Per Update (s)"])
title(["Time by","Number of Agents"])
N = zeros(length(t),nsamps,length(method));
WT = zeros(length(t),nsamps,length(method));
for mi = 1:length(method)
    for si = 1:nsamps
        N(:,si,mi) = sum(out(si,mi).tracked.NT(1:end-1,:),2);
        WT(:,si,mi) = out(si,mi).tracked.wall_time;
    end
    scatter(ax,N(:,:,mi),WT(:,si,mi),5,method_colors(mi,:),'Marker','.')
    ls(mi) = plot(ax,t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end
set(ax,'FontSize',6*res_factor,'XLim',[0,max(N,[],'all')])

% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

f.Position(4) = 0.5*f.Position(4);
f.Children(2).YTick = f.Children(2).YTick([1,end]);
% savefig(sprintf('figs/%s_%s_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_short',name,data_file),'-depsc')

%%
name = "z";
figure("Name","Z");
ax = gobjects(length(method),3);
for ti = 1:3
    for mi = 1:length(method)
        ax(mi,ti) = subplot(length(method),3,r2c(length(method),3,[mi,ti]));
        hold on
        if mi==1
            title(types(ti))
        end
        if ti==1
            ylabel(out(1,mi).method,'FontWeight','bold')
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y = out(si,mi).tracked.mean_z(:,ti);
            plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
        end
        plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
    end
    normalizeYLims(ax)
end
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'FontSize',20)
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
name = "z_mean";
figure("Name","Average Z");
ax_type = gobjects(1,3);
ax_method = gobjects(1,length(method));
for ti = 1:3
    ax_type(ti) = subplot(2,3,ti);
    hold on
    title(types(ti))
    for mi = 1:length(method)
        if ti==1
            method_subplot_width = 3/length(method)-1;
            ax_method(mi) = subplot(2,3,4+(mi-1)*(1+method_subplot_width) + [0,method_subplot_width]);
            title(method(mi))
            hold on
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y(:,si) = out(si,mi).tracked.mean_z(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax_type(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax_type(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
        ls_method(mi,ti) = plot(ax_method(mi),t,ybar,'Color',type_colors(ti,:),'LineWidth',2,'DisplayName',types(ti));
        patch(ax_method(mi),[t;flip(t)],[ybar-s;flip(ybar+s)],type_colors(ti,:),'EdgeColor','none','FaceAlpha',0.6)
   end
   normalizeYLims([ax_type,ax_method])
end
legend(ax_type(end),ls_type(:,end),'Location','best')
legend(ax_method(end),ls_method(end,:),'Location','best')
set([ax_type,ax_method],'XLim',[0 input.simpars.censor_date/(60*24)])
set([ax_type,ax_method],'FontSize',20)
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

fnew = figure;
fnew.Units = "inches";
res_factor = 4;
fnew.Position = [0 0 4 1.5]*res_factor;
for i = 1:3
    ax(i) = subplot(1,3,i);
    copyobj(ax_type(i).Children,ax(i))
    ax(i).XLabel.String = 'Time (d)';
    title(types(i))
    if i==1
        ylabel(["Normalized Distance to","Blood Vessel"])
    end
end
legend(ax(1),flip(ax(1).Children(2:2:end)),'Location','Best')
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
normalizeYLims(fnew)
set(ax,'FontSize',6*res_factor)
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))

% savefig(sprintf('figs/%s_%s_justtype',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_justtype',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_justtype',name,data_file),'-depsc')

fnew.Position(4) = 0.5*fnew.Position(4);
% savefig(sprintf('figs/%s_%s_justtype_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_justtype_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_justtype_short',name,data_file),'-depsc')

%%
name = "phiD";
figure("Name","phi D");
ax = gobjects(length(method),3);
for ti = 1:3
    for mi = 1:length(method)
        ax(mi,ti) = subplot(length(method),3,r2c(length(method),3,[mi,ti]));
        hold on
        if mi==1
            title(types(ti))
        end
        if ti==1
            ylabel(out(1,mi).method,'FontWeight','bold')
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y = out(si,mi).tracked.phiD_mean(:,ti);
            plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
        end
        plot(t,y,'LineWidth',2,'DisplayName',disp_names(mi))
    end
    normalizeYLims(ax)
end
set(ax,'XLim',[0 input.simpars.censor_date/(60*24)])
set(ax,'FontSize',20)
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
name = "phiD_mean";
figure("Name","Mean phi_D");
ax_type = gobjects(1,3);
ax_method = gobjects(1,length(method));
for ti = 1:3
    ax_type(ti) = subplot(2,3,ti);
    hold on
    title(types(ti))
    for mi = 1:length(method)
        if ti==1
            method_subplot_width = 3/length(method)-1;
            ax_method(mi) = subplot(2,3,4+(mi-1)*(1+method_subplot_width) + [0,method_subplot_width]);
            title(method(mi))
            hold on
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y(:,si) = out(si,mi).tracked.phiD_mean(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax_type(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax_type(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
        ls_method(mi,ti) = plot(ax_method(mi),t,ybar,'Color',type_colors(ti,:),'LineWidth',2,'DisplayName',types(ti));
        patch(ax_method(mi),[t;flip(t)],[ybar-s;flip(ybar+s)],type_colors(ti,:),'EdgeColor','none','FaceAlpha',0.6)
    end
    normalizeYLims([ax_type,ax_method])
end
legend(ax_type(end),ls_type(:,end),'Location','best')
legend(ax_method(end),ls_method(end,:),'Location','best')
set([ax_type,ax_method],'XLim',[0 input.simpars.censor_date/(60*24)])
set([ax_type,ax_method],'FontSize',20)
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
name = "codensity_mean";
f=figure("Name","Mean Codensity");
f.Units = "inches";
res_factor = 2;
f.Position = [0 0 4 1.5]*res_factor;
ax_type = gobjects(1,3);
ax_method = gobjects(1,length(method));
for ti = 1:3
    ax_type(ti) = subplot(2,3,ti);
    hold on
    title(types(ti))
    if ti==1
        ylabel('Codensity')
    end
    for mi = 1:length(method)
        if ti==1
            method_subplot_width = 3/length(method)-1;
            ax_method(mi) = subplot(2,3,4+(mi-1)*(1+method_subplot_width) + [0,method_subplot_width]);
            title(disp_names(mi))
            if mi==1
                ylabel('Codensity')
            end
            hold on
        end
        for si = 1:nsamps
            t = out(si,mi).tracked.T/(60*24);
            y(:,si) = out(si,mi).tracked.codensity(:,ti);
        end
        s = std(y,[],2);
        ybar = mean(y,2);
        ls_type(mi,ti) = plot(ax_type(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
        patch(ax_type(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
        ls_method(mi,ti) = plot(ax_method(mi),t,ybar,'Color',type_colors(ti,:),'LineWidth',2,'DisplayName',types(ti));
        patch(ax_method(mi),[t;flip(t)],[ybar-s;flip(ybar+s)],type_colors(ti,:),'EdgeColor','none','FaceAlpha',0.6)
    end


    normalizeYLims([ax_type,ax_method])
end
legend(ax_type(end),ls_type(:,end),'Location','best')
legend(ax_method(end),ls_method(end,:),'Location','best')
set([ax_type,ax_method],'XLim',[0 input.simpars.censor_date/(60*24)])
set([ax_type,ax_method],'FontSize',6*res_factor)
set([ax_type,ax_method],'XTick',[0 input.simpars.censor_date/(60*24)])
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
if isfield(out(1).tracked,'codensity_by_type')
    name = "codensity_mean_by_type";
    f=figure("Name","Mean Codensity");
    f.Units = "inches";
    res_factor = 2;
    f.Position = [0 0 4 2]*res_factor;
    ax_type = gobjects(3,3);
    ls_type = gobjects(length(method),3,3);
    for ti = 1:3
        for oi = 1:3
            ax_type(oi,ti) = subplot(3,3,r2c(3,3,[oi,ti]));
            hold on
            if oi==1
                if ti~=2
                    title(types(ti))
                else
                    title(["Codensity of",types(ti)])
                end
            end
            if ti==1
                if oi~=2
                    ylabel(types(oi),'FontWeight','bold')
                else
                    ylabel(["Relative to",types(oi)],'FontWeight','bold')
                end
            end
            for mi = 1:length(method)
                for si = 1:nsamps
                    t = out(si,mi).tracked.T/(60*24);
                    y(:,si) = out(si,mi).tracked.codensity_by_type(:,oi,ti);
                end
                s = std(y,[],2);
                ybar = mean(y,2);
                ls_type(mi,oi,ti) = plot(ax_type(oi,ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
                patch(ax_type(oi,ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
            end
        end


    end
    normalizeYLims(ax_type)
    legend(ax_type(end),ls_type(:,end,end),'Location','best')
    set(ax_type,'XLim',[0 input.simpars.censor_date/(60*24)])
    set(ax_type,'FontSize',6*res_factor)
    set(ax_type,'XTick',[0 input.simpars.censor_date/(60*24)])
%     savefig(sprintf('figs/%s_%s',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')
end

%%
if isfield(out(1).tracked,'codensity_by_type')
    name = "codensity_of_stem_mean";
    f=figure("Name","Mean Codensity");
    f.Units = "inches";
    res_factor = 2;
    f.Position = [0 0 4 1.5]*res_factor;
    ax_type = gobjects(1,3);
    ls_type = gobjects(length(method),3);
    t = out(1).tracked.T/(60*24);
    y = zeros(length(t),3,nsamps,length(method));
    for ti = 1:3
        ax_type(ti) = subplot(1,3,ti);
        hold on
        title(types(ti))
        if ti==1
            ylabel(["Codensity Relative","to Stem Cells"])
        end
        for mi = 1:length(method)
            for si = 1:nsamps
                y(:,ti,si,mi) = out(si,mi).tracked.codensity_by_type(:,1,ti);
            end
            s = std(y(:,ti,:,mi),[],3);
            ybar = mean(y(:,ti,:,mi),3);
            ls_type(mi,ti) = plot(ax_type(ti),t,ybar,'Color',method_colors(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
            patch(ax_type(ti),[t;flip(t)],[ybar-s;flip(ybar+s)],method_colors(mi,:),'EdgeColor','none','FaceAlpha',0.6)
        end
    end

    normalizeYLims(ax_type)
    legend(ax_type(end),ls_type(:,end,end),'Location','best')
    set(ax_type,'XLim',[0 input.simpars.censor_date/(60*24)])
    set(ax_type,'FontSize',6*res_factor)
%     savefig(sprintf('figs/%s_%s',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

    f=figure;
    f.Units = "inches";
    res_factor = 2;
    f.Position = [0 0 4 1.5]*res_factor;
    ax_method = gobjects(1,length(method));
    ls_method = gobjects(length(method),3);
    for mi = 1:length(method)
        ax_method(mi) = subplot(1,length(method),mi);
        hold on
        for ti = 1:3
            s = std(y(:,ti,:,mi),[],3);
            ybar = mean(y(:,ti,:,mi),3);
            ls_method(mi,ti) = plot(ax_method(mi),t,ybar,'Color',type_colors(ti,:),'LineWidth',2,'DisplayName',types(ti));
        end
        title(disp_names(mi))
        xlabel('Time (d)')
        if mi == 1
            ylabel(["Codensity Relative","to Stem Cells"])
        end
    end
    legend(ax_method(1),ls_method(1,:),'Location','best')
    set(ax_method,'FontSize',6*res_factor)
    set(ax_method,'XLim',[0 input.simpars.censor_date/(60*24)])
    set(ax_method,'XTick',0:10:input.simpars.censor_date/(60*24))

%     savefig(sprintf('figs/%s_%s_bymethod',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s_bymethod',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s_bymethod',name,data_file),'-depsc')

end

%%
name = "il6_by_z";
figure("Name","IL-6 by Z")
ax = gobjects(length(method),nsamps);
for mi = 1:length(method)

    for si = 1:nsamps

        ax(mi,si) = subplot(length(method),nsamps,r2c(length(method),nsamps,[mi,si]));

        contourf(out(si,mi).tracked.T/(60*24),linspace(0,1,out(si,mi).grid.size(3)),out(si,mi).tracked.substrate_by_z(:,:,1)')
        colorbar
    end
end
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

%%
name = "il6_by_z_mean";
f = figure("Name","Average IL-6 by Z");
f.Units = "inches";
res_factor = 2;
f.Position = [0 0 4 1.5]*res_factor;
ax = gobjects(length(method)+1,1);
Z = cell(1,length(method));
ncolors = 15;
n_percent_colors = 100; % number of colors to use in + %s (same for - %s)
for mi = 1:length(method)
    ax(mi) = subplot(1,length(method)+1,mi);
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

cmap = flipud(cubehelix(ncolors,3,-.1,1,1,[0,1],[0.3,1]));
colormap(cmap)

ax(end) = subplot(1,length(method)+1,length(method)+1);
percent_diffs = 100*(mean(Z{length(method)},3)./mean(Z{1},3)-1);
percent_diffs(isnan(percent_diffs)) = 0;
max_percent_diff = max(abs(percent_diffs),[],'all');
[bounds_percent_diff(1),bounds_percent_diff(2)] = bounds(percent_diffs(:));

if isequal(bounds_percent_diff,[0,0])
    bounds_percent_diff = [-1 1];
end
cmap_pd_pos = cubehelix(n_percent_colors,3,-.1,1,1,[0,1],[0.3,1]);
cmap_pd_neg = cubehelix(n_percent_colors,3,+.1,1,1,[0,1],[0.3,1]);

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

set(ax,'FontSize',6*res_factor)
for i = 1:numel(ax)
    ax(i).XLabel.String = 'Time (d)';
end
set(ax,'XTick',0:10:input.simpars.censor_date/(60*24))

% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

f.Position(4) = 0.5*f.Position(4);
delete(f.Children(strcmp(get(f.Children,'type'),'colorbar')))
% savefig(sprintf('figs/%s_%s_short',name,data_file))
% print(sprintf('figs/jpgs/%s_%s_short',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s_short',name,data_file),'-depsc')

%% determine if next two blocks need to be run
ail6r_present = false;
for out_i = 1:numel(out)
    if out(1).solver.is_present(2)
        ail6r_present = true;
        continue;
    end
end

%%
if ail6r_present
    name = "ail6r_by_z";
    figure("Name","aIL-6R by Z")
    ax = gobjects(length(method),nsamps);
    for mi = 1:length(method)

        for si = 1:nsamps

            ax(mi,si) = subplot(length(method),nsamps,r2c(length(method),nsamps,[mi,si]));

            contourf(out(si,mi).tracked.T/(60*24),linspace(0,1,out(si,mi).grid.size(3)),out(si,mi).tracked.substrate_by_z(:,:,2)')
            colorbar
        end
    end
%     savefig(sprintf('figs/%s_%s',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')
end
%%
if ail6r_present
    name = "ail6r_by_z_mean";
    f = figure("Name","Average aIL-6R by Z");
    f.Units = "inches";
    res_factor = 2;
    f.Position = [0 0 4 1.5]*res_factor;
    ax = gobjects(length(method)+1,1);
    Z = cell(1,length(method));
    ncolors = 15;
    n_percent_colors = 100; % number of colors to use in + %s (same for - %s)
    for mi = 1:length(method)
        ax(mi) = subplot(1,length(method)+1,mi);
        for si = 1:nsamps
            Z{mi}(:,:,si) = out(si,mi).tracked.substrate_by_z(:,:,2)';
        end

        contourf(out(si,mi).tracked.T/(60*24),linspace(0,1,out(si,mi).grid.size(3)),mean(Z{mi},3),ncolors,"LineColor","none")
%         LL{mi} = ax(mi).Children.LevelList;
        c=colorbar;
        c.Label.String = 'aIL-6R Concentration (nM)';
        c.Location = "southoutside";
        title(ax(mi),disp_names(mi))
        if mi==1
            ylabel(["Normalized Distance to","Blood Vessel"])
        end
        
    end

    


    cmap = flipud(cubehelix(ncolors,3,-.1,1,1,[0,1],[0.3,1]));
    colormap(cmap)

    ax(end) = subplot(1,length(method)+1,length(method)+1);
    percent_diffs = 100*(mean(Z{length(method)},3)./mean(Z{1},3)-1);
    percent_diffs(isnan(percent_diffs)) = 0;
    max_percent_diff = max(abs(percent_diffs),[],'all');
    [bounds_percent_diff(1),bounds_percent_diff(2)] = bounds(percent_diffs(:));

    if isequal(bounds_percent_diff,[0,0])
        bounds_percent_diff = [-1 1];
    end
    
    cmap_pd_pos = cubehelix(n_percent_colors,3,-.1,1,1,[0,1],[0.3,1]);
    cmap_pd_neg = cubehelix(n_percent_colors,3,+.1,1,1,[0,1],[0.3,1]);

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

%     savefig(sprintf('figs/%s_%s',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')

    f.Position(4) = 0.5*f.Position(4);
    delete(f.Children(strcmp(get(f.Children,'type'),'colorbar')))
%     savefig(sprintf('figs/%s_%s_short',name,data_file))
%     print(sprintf('figs/jpgs/%s_%s_short',name,data_file),'-djpeg')
%     print(sprintf('figs/epss/%s_%s_short',name,data_file),'-depsc')
end