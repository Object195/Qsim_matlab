classdef Q_operator
    %Q_operator A matrix representation of a quantum state or operator
    %  (only square matrices are included)
    
    properties
        matrix = [] %matrix representation 
        dims = {} %dimension of Hilbert space
        qsize %dimension of representation
        %compute the following quantities only if specific functions are
        %called
        Emat = [] %eigenvalue matrix 
        Umat = [] %unitary transformation matrix from eigenbasis to orignal basis
    end
    
    methods
        % Constructor
        function obj = Q_operator(matrix,dims)
            if nargin > 0
                obj.matrix = matrix;
                if nargin >= 2
                    obj.dims = dims;
                %initialize size
                    obj = obj.update_size();
                end
            end
        end
        %update size
        function obj = update_size(obj)
                 fe = cellfun(@(x) x(1), obj.dims); % Extracts the first element of each list
                 obj.qsize = prod(fe); % Multiplies all the first elements together
        end
        % Check compatibility for binary operation
        function isCompatible = check_mat_comp(obj, b)
            %check for two matrix
            isCompatible = all(obj.qsize == b.qsize);
        end

        function isCompatible = check_vec_comp(obj, b)
            %check matrix and vector
            if isa(obj, 'q_rep.Q_operator') && isa(b, 'q_rep.Q_ket')
                isCompatible = all(obj.qsize == b.qsize(1));
            elseif isa(obj, 'q_rep.Q_ket') && isa(b, 'q_rep.Q_operator')
                 isCompatible = all(obj.qsize(2) == b.qsize);
            end
        end
        %copy dimension and size
        function obj = copy_dim(obj,B)
            obj.dims = B.dims;
             obj.qsize = B.qsize;
        end
                
        % Scalar or Matrix addition
        function result = plus(obj, b)
            if isa(obj, 'q_rep.Q_operator') && isa(b, 'q_rep.Q_operator')
                if obj.check_mat_comp(b)
                    result = q_rep.Q_operator(obj.matrix + b.matrix,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_operator') && isnumeric(b)
                result = q_rep.Q_operator(obj.matrix + b,obj.dims);
            elseif isa(b, 'q_rep.Q_operator') && isnumeric(obj)
                result = q_rep.Q_operator(obj + b.matrix,b.dims);
            else
                error('Operand must be a scalar or Q_operator');
            end
        end
        
        % Scalar or Matrix subtraction
        function result = minus(obj, b)
            if isa(obj, 'q_rep.Q_operator') && isa(b, 'q_rep.Q_operator')
                if obj.check_mat_comp(b)
                    result = q_rep.Q_operator(obj.matrix - b.matrix,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_operator') && isnumeric(b)
                result = q_rep.Q_operator(obj.matrix - b,obj.dims);
            elseif isa(b, 'q_rep.Q_operator') && isnumeric(obj)
                result = q_rep.Q_operator(obj-b.matrix,b.dims);
            else
                error('Operand must be a scalar or Q_operator');
            end
        end
        
        % Scalar multiplication or Vector-Matrix multiplication
        function result = mtimes(obj, b)
            if isa(obj, 'q_rep.Q_operator') && isa(b, 'q_rep.Q_operator')
                if obj.check_mat_comp(b)
                    result = q_rep.Q_operator(obj.matrix * b.matrix,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_operator') && isa(b, 'q_rep.Q_ket')
                if obj.check_vec_comp(b)
                    result = q_rep.Q_ket(obj.matrix * b.vec,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_operator') && isnumeric(b)
                result = q_rep.Q_operator(obj.matrix * b,obj.dims);
            elseif isa(b, 'q_rep.Q_operator') && isnumeric(obj)
                result = q_rep.Q_operator(obj*b.matrix,b.dims);
            else
                error('Operand must be a scalar or Q_operator or Q_ket');
            end
        end
        
        % Conjugate transpose
        function result = dag(obj)
            result = q_rep.Q_operator(obj.matrix');
            result = copy_dim(result,obj);
        end
        % Compute Eigenvalue
        function obj = e_decomp(obj)
            if isempty(obj.Emat)
                [obj.Umat,obj.Emat] = eig(obj.matrix);
            end
        end
        % display eigenvalue
        function result = e_val(obj)
            obj = obj.e_decomp();
            result = diag(obj.Emat);
        end
        % exponential of the operator (times some scalar parameter a)
        function result = exp(obj,a)
            r_mat = (obj.Umat)*diag(exp(a*obj.e_val))*ctranspose(obj.Umat);
            result = q_rep.Q_operator(r_mat,obj.dims);
        end
        % Von-Neumann entropy
        function result = entropy(obj)
            result = 0;
            evec = obj.e_val;
            for i = 1:length(evec)
                if evec(i) > 0
                    result = result - evec(i)*log(evec(i));
                end
            end
        end
        %linear entropy
        function result = linear_entropy(obj)
            N = obj.qsize; rho2 = obj*obj;
            result = N/(N-1) * (1-rho2.tr());
        end
        %trace
        function result = tr(obj)
            result = trace(obj.matrix);
        end
        %partial trace
        function result = ptrace(obj,index)
            % index is an array that contains the index of Hilbert spaces
            % to be traced out, returned dm will be the reduced state aranged 
            %with space index in increasing order
            index = sort(index);
            mask = zeros(1,length(obj.dims));
            mask(index) = 1;
            %compute the times of projection needs to be applied
            
            trace_dim = cellfun(@(x) x(1), {obj.dims{index}});
            %generate an index array for constructing all basis
            comb_cell_vec =  arrayfun(@(x) 1:x,trace_dim,'UniformOutput', false);
            %construct all possible indexes for iteration 
            grids = cell(1, length(trace_dim));
            [grids{:}] = ndgrid(comb_cell_vec{:});
            combinations = cellfun(@(x) x(:), grids, 'UniformOutput', false);
            ind_array = [combinations{:}];
            %disp(ind_array)
            %generate identity operator for each subspace (standard basis)
            comb_cell_id =  arrayfun(@(x) eye(x),trace_dim,'UniformOutput', false);
            %display(comb_cell_id)
            %iterate and apply all projections
            result_mat = 0;
            for i = 1:size(ind_array,1) % index for row of iteration grid, index of projection
                %generate transformation matrix
                 j=1; %this index counts subspace to be traced
                 for k = 1:length(obj.dims) %index for all subspace
                     new_d = obj.dims{k}(1);
                    if mask(k)
                    %use basis for the subspace
                        newmat = comb_cell_id{j}(:,ind_array(i,j)); %index for basis of subspace j
                        j = j+1; %move the index j 1 step forward
                    else
                        newmat = eye(new_d);
                    end
                        
                    if k == 1
                        t_mat = newmat;
                    else
                        t_mat = kron(t_mat,newmat);
                    end
                 end
                 result_mat = result_mat + transpose(t_mat)*obj.matrix*t_mat;
            end
            result = q_rep.Q_operator(result_mat,obj.dims(logical(not(mask))));
        end
                       
    end
end

