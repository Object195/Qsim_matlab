a=0.5*[1,0,0,1;0,0,0,0;0,0,0,0;1,0,0,1];
A = q_rep.Q_operator(a,{[2,2],[2,2]});
B=A.ptrace([2]);
B.dims;
%A = q_rep.Q_operator(a,{[2,2]})
%B=q_func.tensor(A,A) 
%B.matrix
v = sqrt(1/3)*func.spin.z_basis(2,[0,0])+sqrt(2/3)*func.spin.z_basis(2,[1,1]);
func.spin.sq_dir(v,1)
            