classdef constant
    % constant matrix and scalars
    
    properties (Constant)
        %Pauli matrices in z representation
        sigma_z = [1,0;0,-1];
        sigma_x = [0,1;1,0];
        sigma_y = [0,1;1,0];
        %basis of spin 1/2
        spin_up = [1;0];
        spin_down = [0;1];
    end
    
   
end

