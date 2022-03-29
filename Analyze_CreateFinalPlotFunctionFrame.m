clearvars;
close all

data_file = "il6_cohort001"; % for example
load(sprintf('data/%s.mat',data_file))

out(1).plot_properties.plotLocations = true;
out(1).plot_properties = initializeFigure_IL6(out(1).plot_properties);
plotFunction_IL6(out(1).plot_properties,40/100,out(1).tumors,out(1).grid,out(1).tracked,length(out(1).tracked.T),out(1).pars.dt,[],out(1).inds)

mf = out(1).plot_properties.my_fig;
%%
f = mf.handle;
f.Units = "inches";
res_factor = 4;
f.Position = [0 0 4 3]*res_factor;
mf.ax(1).Legend.Visible = "off";
for i = 2:5
    mf.ax(i).XAxis.TickValues = 0:10*60*24:input.simpars.censor_date;
    mf.ax(i).XAxis.TickLabels = 0:10:input.simpars.censor_date/(60*24);
    mf.ax(i).XAxis.Limits = [0 input.simpars.censor_date];
    mf.ax(i).FontSize = 6*res_factor;
end
mf.ax(5).XAxis.Label.String = 'Time (d)';
mf.ax(2).XAxis.Label.String = 'Time (d)';
for j = 1:length(mf.ax(1).Children)
    mf.ax(1).Children(j).SizeData = 10*res_factor;
end
axis(mf.ax(1),'square')
view(mf.ax(1),[cos(pi/4+0.23),sin(pi/4+0.23),cos(pi/4)])
mf.ax(2).Legend.String{1}='Stem';
mf.ax(2).Legend.String{2}='Progenitor';

mf.ax(4).Title.String = ["Distance to Blood Vessel"];

y_start = mf.ax(5).Position(2);
scale_factor = (sum(mf.ax(2).Position([2,4]))-mf.ax(5).Position(2))/(sum(mf.ax(2).Position([2,4]))-mf.ax(4).Position(2));
gap = mf.ax(2).Position(2) - sum(mf.ax(3).Position([2,4]));

if ~exist("permuted_panels","var") || ~permuted_panels
    mf.ax(1).Position(4) = mf.ax(1).Position(3);
    mf.ax(1).Position(2) = 1-mf.ax(1).Position(4)-(1-sum(mf.ax(2).Position([2,4])));

    for i = 5:-1:3
        mf.ax(i).Position(2) = mf.ax(i-1).Position(2);
    end
    mf.ax(2).Position(1) = mf.ax(1).Position(1);
    mf.ax(2).Position(2) = mf.ax(5).Position(2);

    for i = 1:5
        mf.ax(i).Position(4) = mf.ax(i).Position(4)*scale_factor;
    end
    permuted_panels = true;
end

%%

mf.ax(1).Position(2) = y_start+mf.ax(2).Position(4)+gap*scale_factor;
mf.ax(2).Position(2) = y_start;
mf.ax(5).Position(2) = y_start;
mf.ax(4).Position(2) = y_start+mf.ax(5).Position(4)+gap*scale_factor;
mf.ax(3).Position(2) = y_start+2*mf.ax(5).Position(4)+2*gap*scale_factor;

name = 'sim_snapshot';
% savefig(sprintf('figs/%s%s_%s',name,out(1).method,data_file))
% print(sprintf('figs/jpgs/%s%s_%s',name,out(1).method,data_file),'-djpeg')
% print(sprintf('figs/epss/%s%s_%s',name,out(1).method,data_file),'-depsc')
