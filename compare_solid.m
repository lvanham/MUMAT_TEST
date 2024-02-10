close all

%% Nearest neighbor effect
ind = 3;
fs = 50;
% tet = '2625';
% p = ['Solid sphere\S' tet '\B_S' tet '_'];
% plot_polerr_solid(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'], [p '75.dat'], [p 'half.dat'])

tet = '5510';
p = ['Solid sphere\S' tet '\B_S' tet '_'];
plot_polerr_solid(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'])
% tet = '10294';
% p = ['Solid sphere\S' tet '\B_S' tet '_'];
% plot_polerr_solid(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'], [p '75.dat'], [p 'half.dat'])


h = findall(groot, 'Type','figure');

for i = 1:length(h)
    CA = h(i).CurrentAxes;
    h(i).OuterPosition = [1 1 1920 1440];
    CA.Legend.FontSize = fs;
    CA.Legend.Location = 'northeast';
    CA.Legend.NumColumns = 1;
    CA.FontSize = fs;
    l = 6;
    for j = 1:l
        h(i).CurrentAxes.Children(j).LineWidth = 12;
    end
   % h(i).CurrentAxes.Children = [h(i).CurrentAxes.Children(2:end);h(i).CurrentAxes.Children(1)];
end

f = gca;
CA.Children(2).DisplayName = 'NN=80%';
CA.Children(3).DisplayName = 'NN=85%';
CA.Children(4).DisplayName = 'NN=90%';
CA.Children(5).DisplayName = 'NN=95%';
CA.Children(6).DisplayName = 'NN=100%';
f.Title.String = 'Solid sphere magnetic field strength error';
f.Title.Interpreter = 'none';
f.YLabel.String = 'Magnetic field strength error';
f.YTick = -0.5:0.1:0.5;
f.Title.FontSize = fs;
