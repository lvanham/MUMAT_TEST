function [] = plot_comps_hollow(varargin)
%PLOT_COMPS_HOLLOW Plots B components from mumat vs polar angle
%   All data must use the grid of points defined in points.dat
%   Call using plot_comps_hollow(ir, path to .dat files )
%   ir: radial index, 1 <= ir <= # of r values
%   path to .dat files: character array, e.g. 'file.dat'

if nargin <= 1
    disp('Not enough arguments')
    return
end

try
    ir = varargin{1};
    if ~isa(ir, 'double')
        disp('Radial index not of type double')
        return
    end
catch
    disp('Could not read radial index')
    return
end

R = 500;
ir = max(ir,1);

thetas = [0:pi/50:pi]; 
ntheta = length(thetas);

xyz = importdata("C:\Users\lucas\Desktop\MUMAT_TEST\Hollow sphere\points_hollow.dat");
dist = R-norm(xyz(1 + ntheta*(ir-1),:));

filec = nargin-1;
hx = gobjects(1,filec);
hz = gobjects(1,filec);
for i = 1:filec
    file = varargin{i+1};
    try
        B_data = importdata(file);
        
    catch
        disp(['Failed to read file ' file])
        return
    end
    B_data_x = B_data(:,1); B_data_z = B_data(:,3);



    if i == 1
        xfig = figure('Name','x component','NumberTitle','off');
    else
        figure(xfig)
    end
    hx(i) = plot(thetas, B_data_x(1 + (ntheta*(ir-1) ):ntheta*ir), ...
        '-o','DisplayName',file);
    hold on

    if i == 1
        zfig = figure('Name','z component','NumberTitle','off');
    else
        figure(zfig)
    end

    hz(i) = plot(thetas, B_data_z(1 + (ntheta*(ir-1) ):ntheta*ir), ...
        '-o','DisplayName',file);
    hold on    

end
    legend(hz,'box','off','interpreter','none')
    xlabel('polar angle [rad]')
    ylabel('$B_z$','interpreter','latex')
    xlim([0 pi])
    %ylim([0 0.05])
    set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
    'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
    'TickLabelInterpreter','latex')
    title(['dx = ' char(string(dist*10)) ' mm'])
    
    figure(xfig)
    legend(hx,'box','off','interpreter','none')
    xlabel('polar angle [rad]')
    ylabel('$B_x$','interpreter','latex')
    xlim([0 pi])
    %ylim([0 1])
    set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
    'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
    'TickLabelInterpreter','latex')
    title(['dx = ' char(string(dist*10)) ' mm'])

end

