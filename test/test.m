%%
%benchmark 2 level system dynamics
H = 2*pi*0.5*func.spin.Pauli_x;
psi0 = func.spin.z_basis(1,[0]);
%compute population of upstate
obs_z = psi0*psi0.dag();
times = 0:0.01:1;
p_z = func.te_solve.se_exact(H,times,psi0,{obs_z});
%% plot result
plot(times, p_z{1});
grid on; % Turn on the grid for better readability

% Label the axes using LaTeX
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14); % Label for the x-axis
ylabel('$p_{z}$', 'Interpreter', 'latex', 'FontSize', 14); % Label for the y-axis




            