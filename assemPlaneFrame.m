function [k_t,node2_used]=assemPlaneFrame(k_t,k_ele,node1,node2, node2_used)
%assemPlaneFrame This function assembles the element stiffness
% matrix k of the plane frame element with nodes
% i and j into the global stiffness matrix K.
% This function returns the global stiffness
% matrix K after the element stiffness matrix
% k is assembled.

global node_type
global di_dof;%节点原始自由度
global se_dof;%累计节点原始自由度


d(1:3)=[se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3];
% if  node_type(node1)==0
%     d(1:3)=se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3;
% elseif node_type(node1)>0
%     d(1:3)=[se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+2,...
%         se_dof(node1)-di_dof(node1)+4+node2_used];
%     node2_used=node2_used+1;
% elseif node_type(node1)<0
%     d(1:3)=[se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+2,...
%         se_dof(node1)-di_dof(node1)+4+node2_used];
%     node2_used=node2_used+1;
% end

if  node_type(node2)==0
    d(4:6)=se_dof(node2)-2:se_dof(node2);%第3类单元有bug
elseif node_type(node2)>0
    d(4:6)=[se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+2,...
        se_dof(node2)-di_dof(node2)+4+node2_used];
    %se_dof(node2)-di_dof(node2)
    node2_used=node2_used+1;
elseif node_type(node2)<0
    d(4:6)=[se_dof(node2)-di_dof(node2)+1:se_dof(node2)-di_dof(node2)+2,...
        se_dof(node2)-di_dof(node2)+4+node2_used];
    %se_dof(node2)-di_dof(node2)
    node2_used=node2_used+1;
end

for ii=1:6
    for jj=1:6
        k_t(d(ii),d(jj))=k_t(d(ii),d(jj))+k_ele(ii,jj);
    end
end
