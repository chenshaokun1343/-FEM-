function [enf,Dnf,node2_used ] = EquivalentNodeForce(L,angle,node1,node2, type, p1, p2, idof, Dnf, node2_used )
% 计算线性分布载荷的等效结点力
% 输入参数:
% L，angle, type ----- 单元信息
% p1 ----- 第一个结点上的分布力集度值
% p2 ----- 第二个结点上的分布力集度值
% idof --- 分布力的种类，它可以是下面几种
% 1 --- 分布轴向力, 以xt为正
% 2 --- 分布横向力, 以yt为正
% 3 --- 分布弯矩, 以逆时针为正
% 返回值:
% enft ----- 局部坐标系下等效结点力向量
% enf ----- 整体坐标系下等效结点力向量\
% dnf ----- 整体坐标系下全结点力向量\

global node_type
global se_dof
global di_dof
enf=zeros(6,1);
enft=zeros(6,1);
switch type
    case 1
    switch idof
        case 1 % 分布轴向力
        enft( 1 ) = (2*p1+p2)*L/6 ;
        enft( 4 ) = (p1+2*p2)*L/6 ;
        case 2 % 分布法向力
        enft( 2 ) = (7*p1+3*p2)*L/20 ;
        enft( 3 ) = (3*p1+2*p2)*L^2/60 ;
        enft( 5 ) = (3*p1+7*p2)*L/20 ;
        enft( 6 ) = -(2*p1+3*p2)*L^2/60 ;
        case 3 % 分布弯矩
        enft( 2 ) = -(p1+p2)/2 ;
        enft( 3 ) = (p1-p2)*L/12 ;%这里我改了符号又改了回来
        enft( 5 ) = (p1+p2)/2 ;
        enft( 6 ) = -(p1-p2)*L/12 ;%这里我改了符号又改了回来
        otherwise
        disp(  '分布力的种类错误，单元方向:%d',L  ) ;
    end
    case 2 
         switch idof
            case 1 % 分布轴向力
            enft( 1 ) = (2*p1+p2)*L/6 ;
            enft( 4 ) = (p1+2*p2)*L/6 ;
            case 2 % 分布法向力
            enft( 2 ) = (16*p1+9*p2)*L/40 ;
            enft( 3 ) = (8*p1+7*p2)*L^2/120 ;
            enft( 5 ) = (4*p1+11*p2)*L/40 ;
            enft( 6 ) = 0;
            enft( 7 ) = 0;
            case 3 % 分布弯矩
            enft( 2 ) = -p1/8-p2*7/8 ;
            enft( 3 ) = -(p1-p2)*L/16 ;
            enft( 5 ) =  p1/8+p2*7/8 ;
            enft( 6 ) = 0 ;
            enft( 7 ) = 0;
            otherwise
            disp(  '分布力的种类错误，单元长:%d',L  ) ;
         end
      case 3 
         switch idof
            case 1 % 分布轴向力
            enft( 1 ) = (2*p1+p2)*L/6 ;
            enft( 4 ) = (p1+2*p2)*L/6 ;
            case 2 % 分布法向力
            enft( 2 ) = (p1+9*p2)*L/2 ;
            enft( 3 ) = (3*p1+5*p2)*L^2/24 ;
            enft( 5 ) = 0 ;
            enft( 6 ) = (p1+3*p2)*L^2/24 ;
            enft( 7 ) = 0;
            case 3 % 分布弯矩
            enft( 2 ) = 0 ;
            enft( 3 ) = p1/3+p2/6 ;
            enft( 5 ) =  0 ;
            enft( 6 ) = p1/6+p2/3 ;
            enft( 7 ) = 0;
            otherwise
            disp(  '分布力的种类错误，单元长:%d',L  ) ;
         end
    case 0
         switch idof
            case 1 % 分布轴向力
            enft( 1 ) = (2*p1+p2)*L/6 ;
            enft( 4 ) = (p1+2*p2)*L/6 ;
            case 2 % 分布法向力
            enft( 2 ) = (2*p1+p2)*L/6 ;
            enft( 3 ) = 0 ;
            enft( 5 ) = (p1+2*p2)*L/6  ;
            enft( 6 ) = 0 ;
            case 3 % 分布弯矩
            enft( 2 ) = (p1+p2)/2 ;
            enft( 3 ) = 0 ;
            enft( 5 ) = 0 ;
            enft( 6 ) = -(p1+p2)/2 ;
            otherwise
            disp(  '分布力的种类错误，单元长:%d',L  ) ;
         end
end
c=cosd(angle);
s=sind(angle);
T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
    0 0 0 -s c 0;
    0 0 0 0 0 1]; %坐标变换

enf = T \ enft(1:6) ; % 把等效结点力转换到整体坐标下
if type==2
    enf(7) = enft(7);
elseif type==3
    enf(7) = enft(7);
end
Dnf(se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3) = ...
        Dnf(se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3) + [enf(1),enf(2),enf(3)]';
if  node_type(node2)==0
    Dnf(se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+3)...
        =Dnf(se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+3) + [enf(4),enf(5),enf(6)]';
elseif node_type(node2)>0
    Dnf(se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+3)...
        =Dnf([se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+2,...
        se_dof(node2)-di_dof(node2)+4+node2_used]) + [enf(4),enf(5),enf(6)]';
    node2_used=node2_used+1;
elseif node_type(node2)<0
    Dnf(se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+3)...
        =Dnf([se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+2,...
        se_dof(node2)-di_dof(node2)+4+node2_used]) + [enf(4),enf(5),enf(6)]';
    node2_used=node2_used+1;
end

if size(enft)==[7,1]
    Dnf(se_dof(node2)-di_dof(node2)+4)...
        =Dnf(se_dof(node2)-di_dof(node2)+4) +enf(7)';
end

return