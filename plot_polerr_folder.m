function [] = plot_polerr_folder(ir, folder)
%PLOT_POLAR_ERROR_FOLDER Plots error in B from mumat vs polar angle
%   All data must use the grid of points defined in points.dat
%   Call using plot_polar_error_folder(ir, folder )
%   ir: radial index, 1 <= ir <= # of r values
%   folder: character array

try
    if ~isa(ir, 'double')
        disp('Radial index not of type double')
        return
    end
catch
    disp('Could not read radial index')
    return
end

ir = max(ir,1);

R = 500;
thetas = [0:pi/50:pi]; 
ntheta = length(thetas);

B_exact = importdata("B_exact.dat");
len = length(B_exact); len = len/ntheta;
if ir > len
    disp(['Radial index out of bounds, clamping to ' num2str(len)])
    ir = len;
end

mod_B_exact = sqrt(sum(B_exact.^2,2));

xyz = importdata("C:\Users\lucas\Desktop\MUMAT_TEST\points.dat");
dist = norm(xyz(1 + ntheta*(ir-1),:))-R;

d = dir(folder);
filec = length(d)-2;
if filec<0
    disp('Folder does not exist.')
    return
elseif filec==0
    disp('Folder is empty.')
    return
end


files = char(d(3:end).name);

h = gobjects(1,filec);
for i = 1:filec

    file = strip(files(i,:));

    try
        B_data = importdata([folder '\' file]);
        mod_B_data = sqrt(sum(B_data.^2,2));
    catch
        disp(['Failed to read file ' folder '\' file])
        return
    end

    if i == 1
        figure('Name','Error with polar angle','NumberTitle','off'); 
    end

    error_modB = (mod_B_data-mod_B_exact)./mod_B_exact;
    
    h(i) = plot(thetas, error_modB(1 + (ntheta*(ir-1) ):ntheta*ir), ...
        '-o','DisplayName',file);
    hold on
    
end
    % plots come out
    yline(0,'k--');
    legend(h,'box','off','interpreter','none')
    xlabel('polar angle [rad]')
    xlim([0 pi])
    ylim([-0.5 0.5])
    ylabel('$(|B|-|B_\mathrm{exact}|)/|B_\mathrm{exact}|$','interpreter','latex')

    set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')

    title(['dx = ' char(string(dist*10)) ' mm'])
end

