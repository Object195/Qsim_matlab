classdef auxi
    %auxiliay functions (non quantum)
    %   此处显示详细说明
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
    end
end

