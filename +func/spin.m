classdef spin
    %Q_FUNCTION 此处显示有关此类的摘要
    %   basic operations of quantum operators and states
    properties (Constant)
        Pauli_x = q_rep.Q_operator(constant.sigma_x,{[2,2]});
        Pauli_y = q_rep.Q_operator(constant.sigma_y,{[2,2]});
        Pauli_z = q_rep.Q_operator(constant.sigma_z,{[2,2]});
        st_0 = q_rep.Q_ket(constant.spin_up,{[2,1]});
        st_1 = q_rep.Q_ket(constant.spin_down,{[2,1]});
        bell_phi_p = 1/sqrt(2)*(func.spin.z_basis(2,[0,0]) ...
            +func.spin.z_basis(2,[1,1]))
        bell_phi_m = 1/sqrt(2)*(func.spin.z_basis(2,[0,0]) ...
            -func.spin.z_basis(2,[1,1]))
        bell_psi_p = 1/sqrt(2)*(func.spin.z_basis(2,[0,1]) ...
            +func.spin.z_basis(2,[1,0]))
        bell_psi_m = 1/sqrt(2)*(func.spin.z_basis(2,[0,1]) ...
            -func.spin.z_basis(2,[1,0]))

    end
    methods (Static)
        function state = z_basis(N,config)
            % create the direct product basis ket for a N spin 1/2 system
            %   Input 
            %   N - int, total number of spin spaces 
            %   config - vec, config of ket, 0 for up and 1 for down
            
            %   Output: Q_ket of direct product
            for j = 1:N
                if config(j) == 0
                    new_s = func.spin.st_0;
                else
                    new_s = func.spin.st_1;
                end
                if j == 1
                    state = new_s;
                else
                    state = func.gen.tensor(state,new_s);
                end
            end
        end

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
                Jop = Jop + func.spin.spin_op(N,j,sub_op);
            end
        end
        function Jop = Jn_op(N,direction)
            %   spin 1/2 operator Jn with specified direction
            %   N - int, total number of spin spaces 
            %   direction- vec, 3D vector for direction
            %   output:
            %    Q_operator
            xop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_x);
            yop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_y);
            zop = func.spin.ammt_tot(N,0.5*func.spin.Pauli_z);
            
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
    end
end

