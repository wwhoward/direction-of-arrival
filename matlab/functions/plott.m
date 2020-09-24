% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [] = plott(est, paths, par, t)
switch par.type
    case '1d'
        figure
        plt = plot(par.phi, est.spectrum);
        hold on; 

        for p=1:length(paths.AoA(:,2)); l_plt(p)=xline(paths.AoA(p,2), '--r'); end
        for e=1:min(length(paths.AoA(:,2)), length(est.peaks_azi)); p_plt(e)=plot(est.peaks_azi(e), est.peak_val(e), 'ro', 'MarkerSize', 10); end
        title(t)
        legend([plt, l_plt(1), p_plt(1)],'Spatial Spectrum', 'True Angle', 'Estimated Angle');
        set(gca,'XTick',0:pi/2:2*pi)
        set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
        xlabel('Azimuth')
    case '2d'
        figure
        plt = surf(par.tht, par.phi, est.spectrum); 
        hold on; 
        shading interp

        for i=1:min(length(paths.AoA(:,2)), length(est.peaks_azi))
            l_plt(i) = plot3(est.peaks_ele(i), est.peaks_azi(i), est.peak_val(i), '.r', 'MarkerSize', 30); 
        end
        for i=1:length(paths.AoA(:,2))
            p_plt(i) = plot3((paths.AoA(i, 1)), paths.AoA(i, 2), max(est.spectrum(:)), '.b', 'MarkerSize',20);
        end
        title(t)
        if ~isempty(est.peaks_azi)
            legend([plt, l_plt(1), p_plt(1)], 'Spatial Spectrum', 'Estimated DoA', 'True DoA');
        else
            legend([plt, p_plt(1)], 'Spatial Spectrum', 'True DoA');
        end
        set(gca,'YTick',0:pi/2:2*pi)
        set(gca, 'XTick', 0:pi/2:pi)
        set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
        set(gca, 'XTickLabel',{'0', 'pi/2','pi'})
        ylabel('Azimuth')
        xlabel('Elevation')
end
end