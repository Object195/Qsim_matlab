classdef Q_space_s
    % Q_SPACE_S this class contains representations of operators in a
    % Hilbert space of N spin 1/2 particles 
    
    
    properties
        N = 1;
        dtype = 'sparse'
        J_tot = {};
    end
    
    methods
        function obj = Q_space_s(N,dtype)
            %Q_SPACE_S constructor
            obj.N = N;
            obj.dtype = dtype;
        end
        
        function obj = update_ammt(obj)
            %    compute all matrix representation of total ammt operators
            %    {Jx,Jy,Jz}
            xop = obj.ammt_tot(0.5*func.spin.Pauli_x);
            yop = obj.ammt_tot(0.5*func.spin.Pauli_y);
            zop = obj.ammt_tot(0.5*func.spin.Pauli_z);
            obj.J_tot = {xop,yop,zop};
        end
        function iop = I(obj)
            %    Identity operator on the spin space
            qsize = 2^obj.N; dims = repmat({[2,2]}, 1, obj.N);
            iop = func.gen.identity(qsize,obj.dtype,dims);
        end
            
        function state = z_basis(obj,config)
            %   create the direct product basis ket for a N spin 1/2 system
            %   Input 
            %   config - vec, config of ket, 0 for up and 1 for down
            
            %   Output: Q_ket of direct product
            
            st_0 = func.spin.st_0.dtype_conv(obj.dtype);
            st_1 = func.spin.st_1.dtype_conv(obj.dtype);
            for j = 1:obj.N
                if config(j) == 0
                    new_s = st_0;
                else
                    new_s = st_1;
                end
                if j == 1
                    state = new_s;
                else
                    state = func.gen.tensor(state,new_s);
                end
            end
        end
        
        function sop = subop(obj,i,sub_op)
            % create the spin operator acting on space i of N spin spaces
            % the specific form is determined by Q_operator sub_op, the
            % dimension of each spin space is determined by sub_op
            %
            %   Input 
            %   i - int, index of spin space to be acting on 
            %   sub_op Q_operator, representation of subspace spin operator
            mask = zeros(obj.N); mask(i) = 1;
            sub_op = sub_op.dtype_conv(obj.dtype);
            for j = 1:obj.N
                if mask(j)
                    new_op = sub_op;
                else
                    new_op = func.gen.identity(sub_op.qsize,obj.dtype);
                end
                if j == 1
                    sop = new_op;
                else
                    sop = func.gen.tensor(sop,new_op);
                end
            end
        end
        
        function Jop = ammt_tot(obj,sub_op)
            %   create the total spin operator with specified dimension
            %   N - int, total number of spin spaces 
            %   sub_op- Q_operator, representation of single subspace spin operator
            %   output:
            %   Q_operator

            for j = 1:obj.N
                if j == 1
                    Jop = obj.subop(j,sub_op);
                else
                    Jop = Jop + obj.subop(j,sub_op);
                end
            end
        end

        function Jop = Jn_op(obj,direction)
            %   spin 1/2 operator Jn with specified direction
            %   direction- vec, 3D vector for direction
            %   Jtot - Cell of Q_operator, {Jx,Jy,Jz}
            %   output:
            %    Q_operator
            
            Jop = (direction(1)*obj.J_tot{1} + direction(2)*obj.J_tot{2} ...
                    +direction(3)*obj.J_tot{3});
        end
        
    end
end

