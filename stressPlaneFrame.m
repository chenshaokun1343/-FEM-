function ElementStress=stressPlaneFrame(A,E,I,L,angle,type,u6,q6)
c=cosd(angle);
s=sind(angle);
C1=A*E/L;
C2=E*I/L^3;
T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
   0 0 0 -s c 0;
    0 0 0 0 0 1]; 
utj=T*u6;
qtj=T*q6;
switch type
    case 1
        ElementStress=[C1*(utj(4)-utj(1))+qtj(1),   6*C2*L*(utj(3)+utj(6))-12*C2*(utj(5)-utj(2))-qtj(2),...
            -4*C2*L*L*utj(3)-2*C2*L*L*utj(6)+6*C2*L*(utj(5)-utj(2))+qtj(3),     C1*(utj(4)-utj(1))+qtj(4),...
            6*C2*L*(utj(3)+utj(6))-12*C2*(utj(5)-utj(2))+qtj(5),    ...
            -4*C2*L*L*utj(6)-2*C2*L*L*utj(3)+6*C2*L*(utj(5)-utj(2))+qtj(6)];
    case 2 
        ElementStress=[C1*(utj(4)-utj(1))+qtj(1),   3*C2*L*(utj(3))-3*C2*(utj(5)-utj(2))-qtj(2),...
            -3*C2*L*L*utj(3)+3*C2*L*(utj(5)-utj(2))+qtj(3),     C1*(utj(4)-utj(1))+qtj(4),...
            3*C2*L*(utj(3))-3*C2*(utj(5)-utj(2))+qtj(5),    qtj(6)];
    case 3 
         ElementStress=[C1*(utj(4)-utj(1))+qtj(1),   -qtj(2),...
            -C2*L*L*utj(3)+qtj(3),     C1*(utj(4)-utj(1))+qtj(4),...
            0,    C2*L*L*utj(3)+qtj(6)];
    case 0
        ElementStress=[C1*(utj(4)-utj(1))+qtj(1),   -qtj(2),...
            0,     qtj(4),...
            qtj(5),   0];

end
