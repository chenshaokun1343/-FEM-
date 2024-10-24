function [Nxdl1,Nxdl2]=ShapeStrainFuncMatrix(xDl,l)
y=0;
for i=1:length(xDl)
    xdl=xDl(i);
    Nxdl1(i,:) = [1-xdl,  0,  0,    xdl,    0,      0];
    Nxdl2(i,:) = [0,  1-3*xdl.^2+2*xdl.^3,  xdl.*l.*(1-xdl).^2,    0,  3*xdl.^2-2*xdl.^3,    xdl.*l.*(xdl-1).*xdl];
    Bxdl1(i,:) = [-1/l,     0,      0,      1/l,        0,      0];
    Bxdl2(i,:) = [0,    6/l^2*(1-2*xdl)*y,  2/l*(2-3*xdl)*y,    0,  6/l^2*(1-2*xdl)*y,  2/l*(2-3*xdl)*y];
    
end