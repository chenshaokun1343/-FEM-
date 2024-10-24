function ksige_ele=PlaneFrameElementGeoStiff(Fp,A,E,I,L,angle,type)
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
C1=A*E/L;
C2=E*I/L^3;

ksige_local=Fp.*[0 , 0 , 0 , 0 , 0 , 0;
    0 , 6/5/L , 1/10 , 0 , -6/5/L , 1/10;
    0 , 1/10 , 2/15/L , 0 , -1/10 , -L/30;
    0 , 0 , 0 , 0 , 0 , 0;
    0 , -6/5/L , *1/10 , 0 , 6/5/L , -1/10;
    0 , 1/10 , -L/30, 0 , -1/10 , 2*L/15];


T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
   0 0 0 -s c 0;
    0 0 0 0 0 1];                %坐标变换矩阵
ksige_ele=T\ksige_local*T;   