classdef boson
    %BOSON operators for bosonic modes
    
    properties
    end
    
    methods (Static)
        function op = up(c)
            % construct creation operator in fock state representation
            % with dimension specified by cutoff c 
            ele = sqrt(1:c)';
            mat = spdiags(ele,-1,c,c);
            op = q_rep.Q_operator(mat,{[c,c]});
        end
        
        function op = down(c)
            % construct anihilation operator in fock state representation
            % with dimension specified by cutoff 
            ele = sqrt(0:c-1)';
            mat = spdiags(ele,1,c,c);
            op = q_rep.Q_operator(mat,{[c,c]});
        end
    end
end

