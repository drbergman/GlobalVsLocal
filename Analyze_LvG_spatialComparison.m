clearvars;
close all

data_file = "il6_cohort001";
load(sprintf('data/%s.mat',data_file))

lm = out(1,1);
gm = out(1,2);
grid = out(1).grid;

inds = input.inds;
events = 0:4;
event_colors = [65,145,79;... % new tumor cells
                121,121,121;... % proliferating cells
                219,59,53;... % apoptotic cells
                95,48,140;... % movement
                0,0,0;... % resting
                68,55,31]/255;
type_markers = {'o','s','^'};
disp_names = ["Stem","Progenitor","TD"];

f = figure;
f.Units = "inches";
res_factor = 4;
f.Position = [0,0,2.5,2.5]*res_factor;

ax(1) = subplot(2,2,1);
tumors = lm.tumors;
num_prolif_attempts = sum(tumors(:,inds.event_ind)==1);
num_prolifs = sum(tumors(:,inds.event_ind)==0);
[event_color_ind,~,~] = find((tumors(:,inds.event_ind)==events)');
prolif_ind = find(event_color_ind==2);
prolif_ind(tumors(prolif_ind,inds.type_ind)==3) = []; % any TD cells that proliferated must have been progenitor cells and so were not contact inhibited
if ~isempty(prolif_ind)
    contact_inhibited_ind = randsample(prolif_ind,num_prolif_attempts-num_prolifs,false);
else
    contact_inhibited_ind = [];
end
event_color_ind(contact_inhibited_ind) = 6;
for ti = 1:3
    I = tumors(:,inds.type_ind)==ti;
    scatter3(tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(1)),...
             tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(2)),...
             tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(3)),...
             10*res_factor,...
             event_colors(event_color_ind(I),:),'filled',...
             'Marker',type_markers{ti});
end
view(ax(1),[cos(pi/4+0.23),sin(pi/4+0.23),cos(pi/4)])
axis(ax(1),[grid.xx([1 end]),grid.yy([1 end]),grid.zz([1 end])])
axis(ax(1),'square')

ax(1).Box = 'off';
ax(1).NextPlot = 'add';
ax(1).Color = ax(1).Parent.Color;
ax(1).XTick = [];
ax(1).YTick = [];
ax(1).ZTick = [];
ax(1).XColor = 'none';
ax(1).YColor = 'none';
ax(1).ZColor = 'none';

title('Local Method','FontSize',6*res_factor);

%%

ax(2) = subplot(2,2,2);
tumors = gm.tumors;
num_prolif_attempts = sum(tumors(:,inds.event_ind)==1);
num_prolifs = sum(tumors(:,inds.event_ind)==0);
[event_color_ind,~,~] = find((tumors(:,inds.event_ind)==events)');
prolif_ind = find(event_color_ind==2);
prolif_ind(tumors(prolif_ind,inds.type_ind)==3) = []; % any TD cells that proliferated must have been progenitor cells and so were not contact inhibited
if ~isempty(prolif_ind)
    contact_inhibited_ind = randsample(prolif_ind,num_prolif_attempts-num_prolifs,false);
else
    contact_inhibited_ind = [];
end
event_color_ind(contact_inhibited_ind) = 6;
for ti = 1:3
    I = tumors(:,inds.type_ind)==ti;
    scatter3(tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(1)),...
             tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(2)),...
             tumors(tumors(:,inds.type_ind)==ti,inds.location_inds(3)),...
             10*res_factor,...
             event_colors(event_color_ind(I),:),'filled',...
             'Marker',type_markers{ti});
end
view(ax(2),[cos(pi/4+0.23),sin(pi/4+0.23),cos(pi/4)])
axis(ax(2),[grid.xx([1 end]),grid.yy([1 end]),grid.zz([1 end])])
axis(ax(2),'square')

ax(2).Box = 'off';
ax(2).NextPlot = 'add';
ax(2).Color = ax(2).Parent.Color;
ax(2).XTick = [];
ax(2).YTick = [];
ax(2).ZTick = [];
ax(2).XColor = 'none';
ax(2).YColor = 'none';
ax(2).ZColor = 'none';

title('Global Method','FontSize',6*res_factor);

%%
type_colors = cubehelix(3,2,-1.5,1,1,[0,.6],[.4,.6]);
% method_style = {'-','-.'};
ax(3) = subplot(2,2,4);
cla;
hold on;
t = out(1).tracked.T/(60*24);
for mi = 1:2
    for si = 1:nsamps
        cod(:,:,mi,si) = out(si,mi).tracked.codensity;
    end
end
cod = mean(cod,4);
mz_percent_diff = 100*(cod(:,:,2)./cod(:,:,1)-1);
for ti = 1:3
    plot(t,mz_percent_diff(:,ti),'LineWidth',0.5*res_factor,'Color',type_colors(ti,:),'DisplayName',disp_names(ti));
end
xlim(ax(3),[0 t(end)])
title(["Percent Difference","of Codensity"])

%%
ax(3).XTick = linspace(0,t(end),3);
xlabel('Time (d)')
ylabel('Percent Difference')
ax(3).FontSize = 6*res_factor;

%%
type_colors = cubehelix(3,2,-1.5,1,1,[0,.6],[.4,.6]);
ax(4) = subplot(2,2,3);
cla;
hold on;
t = out(1).tracked.T/(60*24);
for mi = 1:2
    for si = 1:nsamps
        mz(:,:,mi,si) = out(si,mi).tracked.mean_z;
    end
end
mz = mean(mz,4);
mz_percent_diff = 100*(mz(:,:,2)./mz(:,:,1)-1);
for ti = 1:3
    plot(t,mz_percent_diff(:,ti),'LineWidth',0.5*res_factor,'Color',type_colors(ti,:),'DisplayName',disp_names(ti));
end
xlim(ax(3),[0 t(end)])
title(["Percent Difference","of Distance to BV"])
legend(ax(4),'Location','best')

%%
ax(4).XTick = linspace(0,t(end),3);
xlabel('Time (d)')
ylabel('Percent Difference')
ax(4).FontSize = 6*res_factor;

%%
name = 'sim_snapshot_compare';
% savefig(sprintf('figs/%s_%s',name,data_file))
% print(sprintf('figs/jpgs/%s_%s',name,data_file),'-djpeg')
% print(sprintf('figs/epss/%s_%s',name,data_file),'-depsc')
% print(sprintf('figs/svgs/%s_%s',name,data_file),'-dsvg')
