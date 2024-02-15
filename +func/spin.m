classdef spin
    %Q_FUNCTION 此处显示有关此类的摘要
    %   basic operations of quantum operators and states
    properties (Constant)
        Pauli_x = q_rep.Q_operator(constant.sigma_x,{[2,2]});
        Pauli_y = q_rep.Q_operator(constant.sigma_y,{[2,2]});
        Pauli_z = q_rep.Q_operator(constant.sigma_z,{[2,2]});
        st_0 = q_rep.Q_ket(constant.spin_up,{[2,1]});
        st_1 = q_rep.Q_ket(constant.spin_down,{[2,1]});
    end
    methods (Static)

        function sop = spin_op(N,i,sub_op)
            %create the spin operator acting on space i of N spin spaces
            % the specific form is determined by Q_operator sub_op, the
            % dimension of each spin space is determined by sub_op
            %
            %   Input 
            %   N - int, total number of spin spaces 
            %   i - int, index of spin space to be acting on 
            %   sub_op Q_operator, representation of subspace spin operator
            %   Output: Q_operator of direct product

            mask = zeros(N); mask(i) = 1;
            for j = 1:N
                if mask(j)
                    new_op = sub_op;
                else
                    new_op = func.gen.identity(sub_op.qsize);
                end
                if j == 1
                    sop = new_op;
                else
                    sop = func.gen.tensor(sop,new_op);
                end
            end
        end

        function Jop = ammt_tot(N,sub_op)
            %create the total spin operator with specified dimension
            %   N - int, total number of spin spaces 
            %   sub_op- Q_operator, representation of single subspace spin operator
            %   output:
            %    Q_operator
            Jop = 0;
            for j = 1:N
                Jop = Jop + func.gen.spin_op(N,j,sub_op);
            end
        end
        function vec = MSD_vec(state)
            %   compute the mean spin direction given a state of spin 1/2
            %   systems
            %   state - Q_operator, quantum state of the system 
            %   
            %   output:
            %   3-vec
            N_spin = length(state.dims);
            jx = func.gen.expect(func.spin.ammt_tot(N_spin,func.spin.Pauli_x),state);
            jy = func.gen.expect(func.spin.ammt_tot(N_spin,func.spin.Pauli_y),state);
            jz = func.gen.expect(func.spin.ammt_tot(N_spin,func.spin.Pauli_z),state);
            vec = [jx,jy,jz]/norm([jx,jy,jz]);
        end
    end
end

