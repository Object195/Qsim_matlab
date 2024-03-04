%% basic operations
Hb = q_rep.Q_space_b(2,[10,8],'sparse'); 
ket0 = Hb.fock([0,1]);
adown = Hb.down_m(2);
ket1 = adown*ket0;
rho = ket1*ket1.dag();
rhoA =  rho.ptrace([2]);
rhoB =  rho.ptrace([1]);
disp('state A')
rhoA.matrix
disp('state B')
rhoB.matrix
%%
Ha = q_rep.Q_space_b(1,[5],'sparse');
Ha.up_m(1).matrix
