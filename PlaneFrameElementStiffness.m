function k_ele=PlaneFrameElementStiffness(A,E,I,L,angle,type)
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
% switch type
%     case 1
    %欧拉梁单元刚度矩阵
        k_local=[C1 , 0 , 0 , -C1 , 0 , 0;
                      0 , 12*C2 , 6*C2*L , 0 , -12*C2 , 6*C2*L;
                      0 , 6*C2*L , 4*C2*L*L , 0 , -6*C2*L , 2*C2*L*L;
                      -C1 , 0 , 0 , C1 , 0 , 0;
                      0 , -12*C2 , -6*C2*L , 0 , 12*C2 , -6*C2*L;
                      0 , 6*C2*L , 2*C2*L*L , 0 , -6*C2*L , 4*C2*L*L];    %局部坐标系下的第一类单元刚度矩阵\
%     case 2 
%         k_local=[C1 , 0 , 0 , -C1 , 0 , 0;
%                       0 , 3*C2 , 3*C2*L , 0 , -3*C2 , 0;
%                       0 , 3*C2*L , 3*C2*L*L , 0 , -3*C2*L , 0;
%                       -C1 , 0 , 0 , C1 , 0 , 0;
%                       0 , -3*C2 , -3*C2*L , 0 , 3*C2 , 0;
%                       0 , 0 , 0 , 0 , 0 , 0];               %局部坐标系下的第二类单元刚度矩阵\
%     case 3 
%          k_local=[C1 , 0 , 0 , -C1 , 0 , 0;
%                       0 , 0, 0 , 0 , 0 , 0;
%                       0 , 0 , C2*L*L , 0 , 0, 0;
%                       -C1 , 0 , 0 , C1 , 0 , 0;
%                       0 , 0 , 0 , 0 , 0 , 0;
%                       0 , 0 , 0 , 0 , 0 , -C2*L*L];         %局部坐标系下的第三类单元刚度矩阵\
%     case 0
%         k_local=[C1 , 0 , 0 , -C1 , 0 , 0;
%                       0 , 0, 0 , 0 , 0 , 0;
%                       0 , 0 , 0 , 0 , 0, 0;
%                       -C1 , 0 , 0 , C1 , 0 , 0;
%                       0 , 0 , 0 , 0 , 0 , 0;
%                       0 , 0 , 0 , 0 , 0 , 0];       %局部坐标系下的简支单元刚度矩阵\（可以不使用）
% end
T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
   0 0 0 -s c 0;
    0 0 0 0 0 1];                %坐标变换矩阵
k_ele=T\k_local*T;          %经过坐标变换后的单元刚度矩阵


