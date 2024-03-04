classdef Q_ket
    %Q_ket A vector representation of a quantum state ket/bra
    
    properties
        vec = [] %matrix representation 
        dims = {} %dimension of Hilbert space
        qsize  %dimension of representation
        dtype = '' %data type of matrix representation
    end
    
    methods
        % Constructor
        function obj = Q_ket(vec,dims)
            if nargin > 0
                obj.vec = vec;
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
                 fe1 = cellfun(@(x) x(1), obj.dims); % Extracts the first element of each list
                 fe2 = cellfun(@(x) x(2), obj.dims); % Extracts the second element of each list
                 obj.qsize = [prod(fe1),prod(fe2)]; % Multiplies all the first elements together
        end
        %update data type according to matrix of construction
        function obj = update_dtype(obj)
            if issparse(obj.vec)
                obj.dtype = 'sparse';
            else
                obj.dtype = 'full';
            end
        end
        %change data type
        function obj = dtype_conv(obj,new_type)
            if not(strcmp(new_type, obj.dtype))
                if strcmp(new_type, 'sparse')
                    obj.vec = sparse(obj.vec);
                    obj.dtype = 'sparse';
                elseif strcmp(new_type, 'full')
                    obj.vec = full(obj.vec);
                    obj.dtype = 'full';
                else
                    error('incorrect data type')
                end
            end
        end
        % Addition
        function result = plus(obj, b)
            if isa(obj, 'q_rep.Q_ket') && isa(b, 'q_rep.Q_ket')
                if (obj.qsize == b.qsize)
                    result = q_rep.Q_ket(obj.vec + b.vec,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_ket') && isnumeric(b)
                result = q_rep.Q_ket(obj.vec + b,obj.dims);
            elseif isa(b, 'q_rep.Q_ket') && isnumeric(obj)
                result = q_rep.Q_ket(obj + b.vec,b.dims);
            else
                error('Operand must be a scalar or Q_operator');
            end
        end
        % substraction
        function result = minus(obj, b)
            if isa(obj, 'q_rep.Q_ket') && isa(b, 'q_rep.Q_ket')
                if (obj.qsize == b.qsize)
                    result = q_rep.Q_ket(obj.vec - b.vec,b.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_ket') && isnumeric(b)
                result = q_rep.Q_ket(obj.vec - b,obj.dims);
            elseif isa(b, 'q_rep.Q_ket') && isnumeric(obj)
                result = q_rep.Q_ket(obj - b.vec,b.dims);
            else
                error('Operand must be a scalar or Q_operator');
            end
        end
        %mutliplication
        function result = mtimes(obj, b)
            if isa(obj, 'q_rep.Q_ket') && isa(b, 'q_rep.Q_ket')
                if (obj.qsize(1) == b.qsize(2))&&(obj.qsize(1)>1)
                    new_dims = cellfun(@(x) [x(1),x(1)], obj.dims, 'UniformOutput', false);
                    result = q_rep.Q_operator(obj.vec * b.vec,new_dims);
                elseif (obj.qsize(2) == b.qsize(1))&&(obj.qsize(2)>1)
                    result = obj.vec * b.vec;
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_ket') && isa(b, 'q_rep.Q_operator')
                if (obj.qsize(2) == b.qsize)
                    result = q_rep.Q_ket(obj.vec * b.matrix,obj.dims);
                else
                    error('Dimension does not match');
                end
            elseif isa(obj, 'q_rep.Q_ket') && isnumeric(b)
                result = q_rep.Q_ket(obj.vec * b,obj.dims);
            elseif isa(b, 'q_rep.Q_ket') && isnumeric(obj)
                result = q_rep.Q_ket(obj * b.vec,b.dims);
            else
                error('Operand must be a scalar or Q_operator');
            end
        end
        %transpose conjugate
        function result = dag(obj)
            new_dims = cellfun(@(x) [x(2),x(1)], obj.dims, 'UniformOutput', false);
            result = q_rep.Q_ket(obj.vec',new_dims);
            
        end
    end
 

    
end

