classdef Q_operator
    %Q_operator A matrix representation of a quantum state or operator
    %  (only square matrices are included)
    
    properties
        matrix = [] % sparse matrix representation 
        dims = {} %dimension of Hilbert space
        qsize %dimension of representation
        %compute the following quantities only if specific functions are
        %called
        dtype = '' %data type of matrix representation
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
                    obj = obj.update_dtype();
                end
            end
        end
        %update size
        function obj = update_size(obj)
                 fe = cellfun(@(x) x(1), obj.dims); % Extracts the first element of each list
                 obj.qsize = prod(fe); % Multiplies all the first elements together
        end
        %update data type according to matrix of construction
        function obj = update_dtype(obj)
            if issparse(obj.matrix)
                obj.dtype = 'sparse';
            else
                obj.dtype = 'full';
            end
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
        %change data type
        function obj = dtype_conv(obj,new_type)
            if not(strcmp(new_type, obj.dtype))
                if strcmp(new_type ,'sparse')
                    obj.matrix = sparse(obj.matrix);
                    obj.dtype = 'sparse';
                elseif strcmp(new_type ,'full')
                    obj.matrix = full(obj.matrix);
                    obj.dtype = 'full';
                else
                    error('incorrect data type')
                end
            end
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
        function obj = e_decomp(obj,n)
            %n speicifies the number of eigenvalues/vectors to be computed
            %if n is not provided all eigenvalues will be computed 
            if isempty(obj.Emat)
                if strcmp(n ,'all')
                    n_e = obj.qsize;
                elseif isinteger(n)
                    n_e = n;
                else
                    error('incorrect specification of number of eigenvalues')
                end
                [V,D] = eigs(obj.matrix,n_e);
                if strcmp(obj.dtype,'sparse')
                    obj.Umat = sparse(V);
                    obj.Emat = sparse(1:n_e,1:n_e,diag(D));
                else
                    obj.Umat = V;
                    obj.Emat = D;
                end
            end
        end
        % display eigenvalue
        function result = e_val(obj)
            if strcmp(obj.dtype,'sparse')
                result = spdiags(obj.Emat);
            else
                result = diags(obj.Emat);
            end
        end
        % exponential of the operator (times some scalar parameter a)
        function result = exp(obj,a)
            obj = obj.e_decomp('all');
            if strcmp(obj.dtype,'sparse')
                Dt =  sparse(1:obj.qsize,1:obj.qsize,exp(a*obj.e_val()));
            else
                Dt = diag(exp(a*obj.e_val));
            end
            r_mat = (obj.Umat)*Dt*ctranspose(obj.Umat);
            result = q_rep.Q_operator(r_mat,obj.dims);
        end
        % Von-Neumann entropy
        function result = entropy(obj)
            result = 0;
            obj = obj.e_decomp('all');
            evec = obj.e_val();
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
            %construct all possible indexes for iteration, ind_array
            % In ind_array each row represents the quantum numbers of each subspace 
            % of a basis vector, each column represents a subspace,  
            % example, for a 2 qubit system ind_array reads 
            % [(0,0),(0,1), (1,0),(1,1)]
            
            grids = cell(1, length(trace_dim));
            [grids{:}] = ndgrid(comb_cell_vec{:});
            combinations = cellfun(@(x) x(:), grids, 'UniformOutput', false);
            ind_array = [combinations{:}]; 
            %disp(ind_array)
            %generate identity operator for each subspace (standard basis)
            %iterate and apply all projections
            for i = 1:size(ind_array,1) % row index of iteration grid, index of projection
                %generate transformation matrix
                 j=1; %index for iteration over subspaces to be traced
                 for k = 1:length(obj.dims) %index for iteration over all subspaces
                     new_d = obj.dims{k}(1);
                    if mask(k)
                    %generate std basis (col vector) for the  jth subspace
                    %the non-zero index is given by the (i,j) quantum number 
                        if strcmp(obj.dtype,'sparse')
                            newmat = sparse(ind_array(i,j),1,1,trace_dim(j),1);
                        else
                            newmat = zeros(trace_dim(j),1);
                            newmat(ind_array(i,j)) = 1;
                        end
                        j = j+1; %move the index j 1 step forward
                    else
                        newmat = speye(new_d);
                    end
                        
                    if k == 1
                        t_mat = newmat;
                    else
                        t_mat = kron(t_mat,newmat);
                    end
                 end
                 if i == 1 
                    result_mat = transpose(t_mat)*obj.matrix*t_mat;
                 else
                    result_mat = result_mat + transpose(t_mat)*obj.matrix*t_mat;
                 end
            end
            result = q_rep.Q_operator(result_mat,obj.dims(logical(not(mask))));
        end
                       
    end
end

