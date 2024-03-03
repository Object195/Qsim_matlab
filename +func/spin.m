classdef spin
    %Q_FUNCTION 此处显示有关此类的摘要
    %   basic operations of quantum operators and states
    properties (Constant)
        Pauli_x = q_rep.Q_operator(constant.sigma_x,{[2,2]});
        Pauli_y = q_rep.Q_operator(constant.sigma_y,{[2,2]});
        Pauli_z = q_rep.Q_operator(constant.sigma_z,{[2,2]});
        up = q_rep.Q_operator(constant.up,{[2,2]});
        down = q_rep.Q_operator(constant.down,{[2,2]});
        st_0 = q_rep.Q_ket(constant.spin_up,{[2,1]});
        st_1 = q_rep.Q_ket(constant.spin_down,{[2,1]});
        bell_phi_p = 1/sqrt(2)*(func.spin.z_basis(2,[0,0],'sparse') ...
            +func.spin.z_basis(2,[1,1],'sparse'))
        bell_phi_m = 1/sqrt(2)*(func.spin.z_basis(2,[0,0],'sparse') ...
            -func.spin.z_basis(2,[1,1],'sparse'))
        bell_psi_p = 1/sqrt(2)*(func.spin.z_basis(2,[0,1],'sparse') ...
            +func.spin.z_basis(2,[1,0],'sparse'))
        bell_psi_m = 1/sqrt(2)*(func.spin.z_basis(2,[0,1],'sparse') ...
            -func.spin.z_basis(2,[1,0],'sparse'))

    end
    methods (Static)
        function state = z_basis(N,config,dtype)
            % create the direct product basis ket for a N spin 1/2 system
            %   Input 
            %   N - int, total number of spin spaces 
            %   config - vec, config of ket, 0 for up and 1 for down
            
            %   Output: Q_ket of direct product
           
            st_0 = func.spin.st_0.dtype_conv(dtype);
            st_1 = func.spin.st_1.dtype_conv(dtype);
            for j = 1:N
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

        function sop = spin_op(N,i,sub_op,dtype)
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
            sub_op = sub_op.dtype_conv(dtype);
            for j = 1:N
                if mask(j)
                    new_op = sub_op;
                else
                    new_op = func.gen.identity(sub_op.qsize,dtype);
                end
                if j == 1
                    sop = new_op;
                else
                    sop = func.gen.tensor(sop,new_op,dtype);
                end
            end
        end

        function Jop = ammt_tot(N,sub_op,dtype)
            %create the total spin operator with specified dimension
            %   N - int, total number of spin spaces 
            %   sub_op- Q_operator, representation of single subspace spin operator
            %   output:
            %   Q_operator

            for j = 1:N
                if j == 1
                    Jop = func.spin.spin_op(N,j,sub_op,dtype);
                else
                    Jop = Jop + func.spin.spin_op(N,j,sub_op,dtype);
                end
            end
        end

        function op_cell = J_cell(N,dtype)
            %   generate a cell of Jx, Jy, Jz total ammt operators
            %   N - int, total number of spin spaces 
            %   dtype- str, data type of the matirx representation
            %   output:
            %   cell of Q_operator
            xop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_x,dtype);
            yop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_y,dtype);
            zop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_z,dtype);
            
        end
        
        function Jop = Jn_op(N,direction,dtype)
            %   spin 1/2 operator Jn with specified direction
            %   N - int, total number of spin spaces 
            %   direction- vec, 3D vector for direction
            %   output:
            %    Q_operator

            xop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_x,dtype);
            yop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_y,dtype);
            zop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_z,dtype);
            
            J_mat = (direction(1)*xop.matrix+direction(2)*yop.matrix ...
                    +direction(3)*zop.matrix);
            Jop = q_rep.Q_operator(J_mat,xop.dims);
        end

        

        function vec = MSD_vec(state,normalized)
            %   compute the mean spin direction given a state of spin 1/2
            %   systems
            %   state - Q_operator, quantum state of the system 
            %   
            %   output:
            %   3-vec
            N_spin = length(state.dims);
            jx = func.gen.expect(func.spin.ammt_tot(N_spin,0.5*func.spin.Pauli_x),state);
            jy = func.gen.expect(func.spin.ammt_tot(N_spin,0.5*func.spin.Pauli_y),state);
            jz = func.gen.expect(func.spin.ammt_tot(N_spin,0.5*func.spin.Pauli_z),state);
            if normalized
                vec = [jx,jy,jz]/norm([jx,jy,jz]);
            else
                vec = [jx,jy,jz];
            end
        end
        function result = var_p(state)
            %compute the minimum variance perpendicular to MSD
            %   state - Q_operator/ket, quantum state of the system 
            %   
            %   output:
            %   float
            N_spin = length(state.dims); 
            [theta,phi] = func.auxi.polar_conv(func.spin.MSD_vec(state,1));
            n1 = [-sin(phi),cos(phi),0];
            n2 = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
            J1 = func.spin.Jn_op(N_spin,n1); J2 = func.spin.Jn_op(N_spin,n2);
            J1sq = func.gen.expect(J1*J1,state); 
            J2sq = func.gen.expect(J2*J2,state);
            cov = 0.5*func.gen.expect(J1*J2 + J2*J1, state);
            result  = 0.5*(J1sq+J2sq - sqrt((J1sq-J2sq)^2 + 4*cov^2));       
        end

        function result = sq_dir(state,plot_vec)
            %   compute optimal spin direction
            %   state - Q_operator/ket, quantum state of the system 
            %   plot - bool, if true, plot the vector on Bloch sphere
            %   output:
            %   float
            N_spin = length(state.dims); 
            [theta,phi] = func.auxi.polar_conv(func.spin.MSD_vec(state,1));
            n1 = [-sin(phi),cos(phi),0];
            n2 = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
            J1 = func.spin.Jn_op(N_spin,n1); J2 = func.spin.Jn_op(N_spin,n2);
            J1sq = func.gen.expect(J1*J1,state); 
            J2sq = func.gen.expect(J2*J2,state);
            cov = 0.5*func.gen.expect(J1*J2 + J2*J1, state);
            A = J1sq-J2sq; B = 2*cov;
            nAB = sqrt(A^2+B^2);
            if nAB == 0
                phi_s = 0;
            else
                if B<=0
                    phi_s = 0.5*acos(-A/nAB);
                else
                    phi_s = pi-0.5*acos(-A/nAB);
                end
            end
            result = cos(phi_s)*n1 + sin(phi_s)*n2;
            if plot_vec
                % Create the unit sphere
                [x, y, z] = sphere;
                figure;
                surf(x, y, z, 'FaceAlpha', 0.5); % Set transparency to 0.5
                hold on;
                axis equal; % Ensure the aspect ratio is equal to make the sphere look correct
                grid on;
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                quiver3(0, 0, 0, result(1), result(2), result(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
                title('Optimal squeezing direction');
                hold off;
            end
        end
        function result = sq_w(state)
            %compute the spin squeezing parameter defined by Wineland
            %   state - Q_operator/ket, quantum state of the system 
            %   
            %   output:
            %   float
            N_spin = length(state.dims); 
            var_p = func.spin.var_p(state);
            bnorm = (norm(func.spin.MSD_vec(state,0)))^2;
            result = N_spin*var_p/bnorm;
        end

        function result = sq_ku(state)
            %compute the spin squeezing parameter defined by Kitagawa and
            %Ueda
            %   state - Q_operator/ket, quantum state of the system 
            %   
            %   output:
            %   float
            N_spin = length(state.dims); 
            var_p = func.spin.var_p(state);
            result = 4*var_p/N_spin;
        end
        function result = sq_cell(state_cell,sq_type)
            %compute the spin squeezing parameter for a cell
            %of states of the same dimension
            %   state_cell - cell of Q_operator/ket, quantum state of the system 
            %   sq_type - str, can be 'w' (Wineland) or 'u' (Kitagawa and Ueda) 
            %   output:
            %   1d vector of float
            N = length(state_cell{1}.dims); 
            Jx = func.spin.ammt_tot(N,0.5*func.spin.Pauli_x);
            Jy = func.spin.ammt_tot(N,0.5*func.spin.Pauli_y);
            Jz = func.spin.ammt_tot(N,0.5*func.spin.Pauli_z);
            result = zeros(1,length(state_cell));
            for i = 1:length(state_cell)
                state = state_cell{i};
                %compute MSD
                jx = func.gen.expect(Jx,state);
                jy = func.gen.expect(Jy,state);
                jz = func.gen.expect(Jz,state);
                J_abs = jx^2 + jy^2 + jz^2;
                %compute MSD and perpendicular directions
                [theta,phi] = func.auxi.polar_conv([jx,jy,jz]/sqrt(J_abs));
                n1 = [-sin(phi),cos(phi),0];
                n2 = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
                J1 = Jx*n1(1) +   Jy*n1(2) + Jz*n1(3);
                J2 = Jx*n2(1) +   Jy*n2(2) + Jz*n2(3);
                %compute expectation
                J1sq = func.gen.expect(J1*J1,state); 
                J2sq = func.gen.expect(J2*J2,state);
                cov = 0.5*func.gen.expect(J1*J2 + J2*J1, state);
                J_var  = 0.5*(J1sq+J2sq - sqrt((J1sq-J2sq)^2 + 4*cov^2));
                if sq_type == 'w'
                    result(i) = N * J_var / J_abs;
                elseif sq_type == 'u'
                    result(i) = 4*J_var/N;
                else
                    error('Undefined type of spin-squeezing parameter.')
                end
            end
        end
        function result = sq_z(state_cell,sq_type)
            %compute the spin squeezing parameter for a cell
            %of states that polarized along z direction
            %   state_cell - cell of Q_operator/ket, quantum state of the system 
            %   sq_type - str, can be 'w' (Wineland) or 'u' (Kitagawa and Ueda) 
            %   output:
            %   1d vector of float
            N = length(state_cell{1}.dims); 
            J1 = func.spin.ammt_tot(N,0.5*func.spin.Pauli_x);
            J2 = func.spin.ammt_tot(N,0.5*func.spin.Pauli_y);
            J3 = func.spin.ammt_tot(N,0.5*func.spin.Pauli_z);
            J1sq_op = J1*J1; J2sq_op = J2*J2; 
            cov_op = 0.5 * (J1*J2 + J2*J1);
            result = zeros(1,length(state_cell));
            for i = 1:length(state_cell)
                state = state_cell{i};
                J1sq = func.gen.expect(J1sq_op,state); 
                J2sq = func.gen.expect(J2sq_op,state);
                cov = func.gen.expect(cov_op, state);
                J_abs = (func.gen.expect(J3, state))^2;
                J_var  = 0.5*(J1sq+J2sq - sqrt((J1sq-J2sq)^2 + 4*cov^2));
                if sq_type == 'w'
                    result(i) = N * J_var / J_abs;
                elseif sq_type == 'u'
                    result(i) = 4*J_var/N;
                else
                    error('Undefined type of spin-squeezing parameter.')
                end
            end
        end
    end
end

