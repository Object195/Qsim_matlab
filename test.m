a=0.5*[1,0,0,1;0,0,0,0;0,0,0,0;1,0,0,1];
A = q_rep.Q_operator(a,{[2,2],[2,2]});
B=A.ptrace([2]);
B.dims;
%A = q_rep.Q_operator(a,{[2,2]})
%B=q_func.tensor(A,A) 
%B.matrix
C=func.spin.Pauli_z;
state = func.spin.st_1;
func.gen.var(C,state)
