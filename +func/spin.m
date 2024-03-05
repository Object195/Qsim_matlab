classdef spin
    %Q_FUNCTION 
    %   basic operations of quantum operators and states
    properties (Constant)
        Pauli_x = q_rep.Q_operator(constant.sigma_x,{[2,2]});
        Pauli_y = q_rep.Q_operator(constant.sigma_y,{[2,2]});
        Pauli_z = q_rep.Q_operator(constant.sigma_z,{[2,2]});
        up = q_rep.Q_operator(constant.up,{[2,2]});
        down = q_rep.Q_operator(constant.down,{[2,2]});
        st_0 = q_rep.Q_ket(constant.spin_up,{[2,1]});
        st_1 = q_rep.Q_ket(constant.spin_down,{[2,1]});
        
        %bell states
        
    end
    methods (Static)
        function state = bell(stype)
            %generate bell state
            S_2 = q_rep.Q_space_s(2,'sparse');
            if strcmp(stype, 'phi_p')
                state = 1/sqrt(2)*(S_2.z_basis([0,0])+S_2.z_basis([1,1]));
            elseif strcmp(stype, 'phi_m')
                state = 1/sqrt(2)*(S_2.z_basis([0,0])-S_2.z_basis([1,1]));
            elseif strcmp(stype, 'psi_p')
                state = 1/sqrt(2)*(S_2.z_basis([0,1])+S_2.z_basis([1,0]));
            elseif strcmp(stype, 'psi_m ')
                state = 1/sqrt(2)*(S_2.z_basis([0,1])-S_2.z_basis([1,0]));
            end
        end
        
        function vec = MSD_vec(state,s_space,normalized)
            %   compute the mean spin direction given a state of spin 1/2
            %   systems
            %   state - Q_operator, quantum state of the system 
            %   s_space - Q_space_s object
            %   Jtot - Cell of Q_operator, {Jx,Jy,Jz}
            %   output:
            %   3-vec

            jx = func.gen.expect(s_space.J_tot{1},state);
            jy = func.gen.expect(s_space.J_tot{2},state);
            jz = func.gen.expect(s_space.J_tot{3},state);
            if normalized
                vec = [jx,jy,jz]/norm([jx,jy,jz]);
            else
                vec = [jx,jy,jz];
            end
        end
        function result = var_p(state,s_space)
            %compute the minimum variance perpendicular to MSD
            %   state - Q_operator/ket, quantum state of the system 
            %   Jtot - Cell of Q_operator, {Jx,Jy,Jz}
            %   output:
            %   float
            [theta,phi] = func.auxi.polar_conv(func.spin.MSD_vec(state,s_space,1));
            n1 = [-sin(phi),cos(phi),0];
            n2 = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
            J1 = s_space.Jn_op(n1); J2 = s_space.Jn_op(n2);
            J1sq = func.gen.expect(J1*J1,state); 
            J2sq = func.gen.expect(J2*J2,state);
            cov = 0.5*func.gen.expect(J1*J2 + J2*J1, state);
            result  = 0.5*(J1sq+J2sq - sqrt((J1sq-J2sq)^2 + 4*cov^2));       
        end

        function result = sq_dir(state,s_space,plot_vec)
            %   compute optimal spin direction
            %   state - Q_operator/ket, quantum state of the system 
            %   plot - bool, if true, plot the vector on Bloch sphere
            %   output:
            %   float
            [theta,phi] = func.auxi.polar_conv(func.spin.MSD_vec(state,s_space,1));
            n1 = [-sin(phi),cos(phi),0];
            n2 = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
            J1 = s_space.Jn_op(n1); J2 = s_space.Jn_op(n2);
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
        function result = sq_w(state,s_space)
            %compute the spin squeezing parameter defined by Wineland
            %   state - Q_operator/ket, quantum state of the system 
            %   s_space - Q_space_s object, spin Hilbert space
            %   output:
            %   float
            N_spin = length(state.dims); 
            var_p = func.spin.var_p(state,s_space);
            bnorm = (norm(func.spin.MSD_vec(state,s_space,0)))^2;
            result = N_spin*var_p/bnorm;
        end

        function result = sq_ku(state,s_space)
            %compute the spin squeezing parameter defined by Kitagawa and
            %Ueda
            %   state - Q_operator/ket, quantum state of the system 
            %   
            %   output:
            %   float
            N_spin = length(state.dims); 
            var_p = func.spin.var_p(state,s_space);
            result = 4*var_p/N_spin;
        end

        function result = sq_z(state_cell,s_space,sq_type)
            %compute the spin squeezing parameter for a cell
            %of states that polarized along z direction
            %   state_cell - cell of Q_operator/ket, quantum state of the system 
            %   sq_type - str, can be 'w' (Wineland) or 'u' (Kitagawa and Ueda) 
            %   output:
            %   1d vector of float
            N = length(state_cell{1}.dims); 
            J1 = s_space.J_tot{1}; 
            J2 = s_space.J_tot{2}; 
            J3 =s_space.J_tot{3};
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

