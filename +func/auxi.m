classdef auxi
    %auxiliay functions (non quantum)
    methods(Static)
        function [theta,phi] = polar_conv(vec)
            %   convert a unit vector to polar coordinates
            theta = acos(vec(3));
            if (vec(1)==0)&&(vec(2)==0)
                phi = 0;
            else
                xy_proj = sqrt(vec(1)^2+vec(2)^2);
                if vec(2)>0
                    phi = acos(vec(1)/xy_proj);
                else
                    phi = 2*pi - acos(vec(1)/xy_proj);
                end
            end
        end

        function [s_mean, sq_ratio] = sq_statistic(s_vec, sq_vec, N_grid, plot_result)
            % discretize entanglement entropy and corresponding squeezing into N grids 
            % and compute std/mean of squeezing each grid
            % s_vec, -vector of floats, array of entanglement parameter
            % sq_vec, -vector of floats, array of squeezing parameter
            % N_grid, -int, number of grids
            edges = linspace(0, max(s_vec), N_grid+1);
            [~,~,bin] = histcounts(s_vec, edges);
            s_mean = zeros(1, N_grid);
            sq_ratio = zeros(1, N_grid);
            for i = 1:N_grid
                s_bin = s_vec(bin == i);
                sq_bin = sq_vec(bin == i);
                s_mean(i) = mean(s_bin);
                sq_ratio(i) = std(sq_bin)/mean(sq_bin);
            end
            if plot_result
                figure1 = figure;
                set(figure1, 'Position', [100, 100, 600, 400])
                plot(s_mean,sq_ratio,'.','MarkerSize',25);
                hold on
                plot(linspace(0, max(s_vec), N_grid), 0.1*ones(size(s_mean)), ...
                '--', 'DisplayName', '$\xi = 1$','LineWidth', 1.5,'Color',[0 0 0]);
  
                set(gca,'FontSize',14)
                xlabel('$\bar{S}$','Interpreter','latex','FontSize',20)
                ylabel('$\Delta \xi/\bar{\xi}$','Interpreter','latex','FontSize',20)
                ylim([0,1])
                grid on

            end
        end
    end
end



           

