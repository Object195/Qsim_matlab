classdef gen
    %general quantum operations
    methods(Static)
        function result = expect(obs,state)
            %   compute first moment of observable given a quantum state
            %   Input 
            %   obs - Q_operator, observable to be measured
            %   state - Q_operator/Q_ket, quantum state (ket or density matrix)
            %   Output: real float
            if isa(obs, 'q_rep.Q_operator') && isa(state, 'q_rep.Q_operator')
                C = obs*state; result = trace(C.matrix);
            elseif (isa(obs, 'q_rep.Q_operator') && isa(state, 'q_rep.Q_ket'))
                C = obs*(state*state.dag()); result = trace(C.matrix);
            else
                error('Inputs must be instances of Q_operator and Q_operator/ket.');
            end
        end

        function result = var(obs,state)
            %   compute first moment of observable given a quantum state
             %   Input 
            %   obs - Q_operator, observable to be measured
            %   state - Q_operator/Q_ket, quantum state (ket or density matrix)
            %   Output: real float
            m1 = func.gen.expect(obs,state);
            if isa(obs, 'q_rep.Q_operator') && isa(state, 'q_rep.Q_operator')
                C = obs*obs*state; result = trace(C.matrix)-m1^2;
            elseif (isa(obs, 'q_rep.Q_operator') && isa(state, 'q_rep.Q_ket'))
                C = obs*obs*(state*state.dag()); result = trace(C.matrix)-m1^2;
            else
                error('Inputs must be instances of Q_operator and Q_operator/ket.');
            end
        end
        function iop = identity(N,dims)
            %   create an N dimensional identity operator
            %   Input 
            %   N - int, dimensnion of the space 
            %   
            %   Output: Q_operator 
            if nargin > 1
                iop = q_rep.Q_operator(eye(N),dims);
            else
                iop = q_rep.Q_operator(eye(N),{[N,N]});
            end
        end

        function t_prod = tensor(A,B)
            %   compute the tensor product of two operators/ states
            %   Input 
            %   A,B - Q_operator/Q_ket, two quantum operators/ states
            %   
            %   Output: Q_operator/Q_ket 
            if isa(A, 'q_rep.Q_operator') && isa(B, 'q_rep.Q_operator')
                result_m = kron(A.matrix,B.matrix);
                new_dim = [A.dims,B.dims];
                t_prod = q_rep.Q_operator(result_m,new_dim);
            elseif (isa(A, 'q_rep.Q_ket') && isa(B, 'q_rep.Q_ket'))
                result_m = kron(A.vec,B.vec);
                new_dim = [A.dims,B.dims];
                t_prod = q_rep.Q_ket(result_m,new_dim);
            else
                error('Both inputs must be instances of Q_operator or Q_ket.');
            end
        end


    
    
    
    end
end

