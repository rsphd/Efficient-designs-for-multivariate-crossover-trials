syms r
assume(r,"real");
assumeAlso(r>0);
assumeAlso(r<1);
V1=[1 r r^2; r 1 r; r^2 r 1];
VR=[1 r^2 r^4; r^2 1 r^2; r^4 r^2 1];
J=[1 1 1; 1 1 1; 1 1 1];
one=[1 1 1];
I=[1 0 0; 0 1 0; 0 0 1];
H=I-J/3;
psi=[0 0 0; 1 0 0; 0 1 0];

V1_star=inv(V1) - inv(V1)*J*inv(V1)/(one*inv(V1)*transpose(one));
simplify(V1_star)
simplify(trace(V1_star*psi))
simplify(trace(H*transpose(psi)*V1_star*psi))

VR_star=inv(VR) - inv(VR)*J*inv(VR)/(one*inv(VR)*transpose(one));
simplify(VR_star)
simplify(trace(VR_star*psi))
simplify(trace(H*transpose(psi)*VR_star*psi))