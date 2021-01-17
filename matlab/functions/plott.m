% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [] = plott(estimate, paths, par, stats)
for estimator = par.estimator
    for smoothing = par.smooth
        est = estimate.(estimator).(smoothing);
        stat = stats.(estimator).(smoothing);
        switch par.type
            case '1d'
                figure
                plt = plot(par.phi, est.spectrum);
                hold on;
                
                for p=1:length(paths.AoA(:,2)); l_plt(p)=xline(paths.AoA(p,2), '--r'); end
                for e=1:min(length(paths.AoA(:,2)), length(est.peaks_azi)); p_plt(e)=plot(est.peaks_azi(e), est.peak_val(e), 'ro', 'MarkerSize', 10); end
                title(est.smoothing + " - "+string(par.K)+" Total Paths - " + est.label)
                legend([plt, l_plt(1), p_plt(1)],'Spatial Spectrum', 'True Angle', 'Estimated Angle');
                set(gca,'XTick',0:pi/2:2*pi)
                set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
                xlabel('Azimuth')
            case '2d'
                figure
                plt = surf(par.tht, par.phi, est.spectrum);
                hold on;
                shading interp
                
                if ~isempty(stat.peak_val)
                    for i=1:par.K
                        l_plt(i) = plot3(stat.uh(2,i), stat.uh(1,i), stat.peak_val(max(1, stat.peak_idx(i))), '.r', 'MarkerSize', 30);
                        %l_plt(i) = plot3(est.peaks_ele(i), est.peaks_azi(i), est.peak_val(i), '.r', 'MarkerSize', 30);
                    end
                end
                for i=1:length(paths.AoA(:,2))
                    p_plt(i) = plot3((paths.AoA(i, 1)), paths.AoA(i, 2), max(est.spectrum(:)), '.b', 'MarkerSize',20);
                end
                title(est.smoothing + " - "+string(par.K)+" Total Paths - " + est.label)
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
end
end