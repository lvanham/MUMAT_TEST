close all
fs = 50;
%% Nearest neighbor effect
ind = 981;
% tet = '2516';
% p = ['Hollow sphere\H' tet '\B_H' tet '_'];
% plot_polerr_hollow(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'], [p '75.dat'], [p 'half.dat'])
% 
% tet = '5292';
% p = ['Hollow sphere\H' tet '\B_H' tet '_'];
% plot_polerr_hollow(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'], [p '75.dat'], [p 'half.dat'])

tet = '10085';
p = ['Hollow sphere\H' tet '\B_H' tet '_'];
    l = 5;

plot_polerr_hollow(ind, [p 'full.dat'], [p '95.dat'], [p '90.dat'], [p '85.dat'],[p '80.dat'])


h = findall(groot, 'Type','figure');

for i = 1:length(h)
    CA = h(i).CurrentAxes;

    h(i).OuterPosition = [1 1 1920 1440];
    CA.Legend.FontSize = fs;
    CA.Legend.Location = 'northeast';
    %CA.Legend.NumColumns = l;
    CA.YLim = [-0.5 0.5];
    CA.FontSize = fs;
    for j = 1:length(CA.Children)
        CA.Children(j).LineWidth = 12;
        if j>1
            CA.Children(j).MarkerSize = 18;
        end
    end
    
    %h(i).CurrentAxes.Children = [h(i).CurrentAxes.Children(2:end);h(i).CurrentAxes.Children(1)];
    
end

%%
f = gca;
CA.Children(2).DisplayName = 'NN=80%';
CA.Children(3).DisplayName = 'NN=85%';
CA.Children(4).DisplayName = 'NN=90%';
CA.Children(5).DisplayName = 'NN=95%';
CA.Children(6).DisplayName = 'NN=100%';
f.Title.String = 'Hollow sphere magnetic field strength error';
f.YLabel.String = 'Magnetic field strength error';
f.YTick = -0.5:0.1:0.5;
%f.Position=1e3*[-0.849400000000000   0.978600000000000   1.216800000000000   0.826400000000000];
