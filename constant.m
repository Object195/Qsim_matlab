classdef constant
    % constant matrix and scalars
    
    properties (Constant)
        %Pauli matrices in z representation
        sigma_z = sparse([1,0;0,-1]);
        sigma_x = sparse([0,1;1,0]);
        sigma_y = sparse([0,-1i;1i,0]);
        up = sparse([0,1;0,0])
        down = sparse([0,0;1,0])
        %basis of spin 1/2
        spin_up = sparse([1;0]);
        spin_down = sparse([0;1]);
    end
    
   
end

