classdef te_solve
    %solve time evolution
    methods(Static)

        function result = evol_dm(H,t,rho0)
            %subfunction, compute the evolved density at some time t given 
            %H,t by diagnolization
            U = H.exp(-1i*t); 
            result = U*rho0*U.dag();
        end
            
        function result = se_exact(H,tsample,i_state,obs)
            % Solve the dynamics of a time independent Hamiltonian
            % by computing the exact propagator
            % Input:
            % H- Q_operator, the Hamiltonain of the system 
            % t_sample, -vector, time points for sampling
            % obs, cell of Q_operator, observable 
            H = H.e_decomp('all'); %Diagonalize
            %compute propagator
            if isa(i_state, 'q_rep.Q_ket')
                if isempty(obs)
                    result = arrayfun(@(U) H.exp(-1i*t)*i_state,tsample, ...
                        'UniformOutput', false);
                else
                    result = {};
                    for A = obs
                        A_ex = arrayfun(@(t) func.gen.expect(A{1},H.exp(-1i*t)*i_state),tsample);
                        result = [result,A_ex];
                    end
                end
            elseif isa(i_state, 'q_rep.Q_operator')
                if isempty(obs)
                    %result = arrayfun(@(t) H.exp(-1i*t)*i_state*H.exp(1i*t),tsample, ...
                        %'UniformOutput', false);
                    result = arrayfun(@(t) func.te_solve.evol_dm(H,t,i_state),tsample, ...
                        'UniformOutput', false);
                else
                    result = {};
                    for A = obs
                        A_ex = arrayfun(@(t) func.gen.expect(A{1},func.te_solve.evol_dm(H,t,i_state),tsample), ...
                            tsample);
                        result = [result,A_ex];
                    end
                end
            else
                error('invalid initial state')
            end


        end
    end
end

