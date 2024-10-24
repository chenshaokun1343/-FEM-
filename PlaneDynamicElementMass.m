function m_ele=PlaneDynamicElementMass(A,rou,L,angle)
%PlaneFrameElementStiffness This function returns the element
% stiffness matrix for a plane frame
% element with modulus of elasticity E,
% cross-sectional area A, moment of
% inertia I, length L, and angle
% theta (in degrees).
% The size of the element stiffness
% matrix is 6 x 6.

c=cosd(angle);
s=sind(angle);

%m_local=rou*A*L/2*diag([1,1,1e-10,1,1,1e-10]);
m_local=rou*A*L/420.*[140, 0, 0, 70, 0, 0;
                                    0, 156, 22*L, 0, 54, -13*L;
                                    0, 22*L, 4*L.^2, 0, 13*L, -3*L.^2;
                                    70, 0, 0, 140, 0, 0;
                                    0, 54, 13*L, 0, 156, -22*L;
                                    0, -13*L, -3*L.^2, 0, -22*L, 4*L.^2];

T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
   0 0 0 -s c 0;
    0 0 0 0 0 1];                %坐标变换矩阵
m_ele=T\m_local*T;          %经过坐标变换后的单元刚度矩阵


