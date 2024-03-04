classdef Q_space_b
    %Q_SPACE_B this class contains representations of operators in a
    % Hilbert space of N Bosonic modes
    
    properties
        N = 1; %number of bosonic spaces
        dtype = 'sparse' %data storage type
        cutoff = []; % vector of int, cutoff for each bosonic subspace
    end
    
    methods
        function obj = Q_space_b(N,cutoff,dtype)
            %Q_SPACE_B constructor
            obj.N = N;
            obj.cutoff = cutoff;
            obj.dtype = dtype;
            
        end

        function iop = I(obj)
            %    Identity operator on the bosonic space
            qsize = prod(obj.cutoff); 
            dims = arrayfun(@(x) repmat(x, 1, 2), obj.cutoff, 'UniformOutput', false);
            iop = func.gen.identity(qsize,obj.dtype,dims);
        end

        function state = fock(obj,config)
            % construct a fock basis ket, the excitation number is given in
            % config, which must have the same size as cutoff
            % the excitation number can not exceeds cutoff
            if length(config) == length(obj.cutoff)
                for j = 1:obj.N
                    if strcmp(obj.dtype, 'full')
                        bvec = zeros(obj.cutoff(j),1);
                        bvec(config(j)+1) = 1;
                        new_s = q_rep.Q_ket(bvec,{[obj.cutoff(j),1]});
                    elseif strcmp(obj.dtype, 'sparse')
                        bvec = sparse(config(j)+1, 1, 1, obj.cutoff(j), 1);
                        new_s = q_rep.Q_ket(bvec,{[obj.cutoff(j),1]});
                    end
                    if j == 1
                        state = new_s;
                    else
                        state = func.gen.tensor(state,new_s);
                    end
                end
            else
                error('incorrect config of the state')
            end
        end
        
        function bop = up_m(obj,m)
            % construct the creation operator on the mth mode
            %
            mask = zeros(obj.N); mask(m) = 1;
            sub_op = func.boson.up(obj.cutoff(m)).dtype_conv(obj.dtype);
            for j = 1:obj.N
                if mask(j)
                    new_op = sub_op;
                else
                    new_op = func.gen.identity(obj.cutoff(j),obj.dtype);
                end
                if j == 1
                    bop = new_op;
                else
                    bop = func.gen.tensor(bop,new_op);
                end
            end
        end

        function bop = down_m(obj,m)
            % construct the annihilation operator on the mth mode
            %
            mask = zeros(obj.N); mask(m) = 1;
            sub_op = func.boson.down(obj.cutoff(m)).dtype_conv(obj.dtype);
            for j = 1:obj.N
                if mask(j)
                    new_op = sub_op;
                else
                    new_op = func.gen.identity(obj.cutoff(j),obj.dtype);
                end
                if j == 1
                    bop = new_op;
                else
                    bop = func.gen.tensor(bop,new_op);
                end
            end
        end
    end
end

