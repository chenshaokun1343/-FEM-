clc,clear,close all
node=[-60 0 0;
    -60 120 0;
    60 120 0;
    60 0 0];
% node=[-60 0 0;
%       -60 40 0;
%       -60 80 0
%       -60 120 0;
%       -20 120 0;
%       20 120 0;
%       60 120 0;
%       60 80 0;
%       60 40 0;
%       60 0 0];   %节点位置，行号为节点编号，1~3列分别为x,y,z方向坐标
ele=[1 2 1;
    2 3 1;
    3 4 1];
% ele=[1 2 1;
%     2 3 1;
%     3 4 1;
%     4 5 1;
%     5 6 1;
%     6 7 1;
%     7 8 1;
%     8 9 1;
%     9 10 1];        %单元信息，行号为单元编号，1~2各列为单元上的节点号码,3列为基本构型
%记得改f(10)=0;
E=3e7.*ones(1,3);             %弹性模量
I=[200.*ones(1,1) 100.*ones(1,1) 200.*ones(1,1)];           %杆截惯性矩
A=[2.*ones(1,1)  1.6.*ones(1,1)  2.*ones(1,1)];              %杆截面面积
rou=1e3.*ones(1,3);      %杆密度
n_node=length(node(:,1));%节点数
n_ele=length(ele(:,1));   %单元数
global node_type
node_type=zeros(1,n_node);
node_type( ele(find(ele(:,3)==2),2) )=length(ele(find(ele(:,3)==2)));
node_type( ele(find(ele(:,3)==3),2) )=-length(ele(find(ele(:,3)==3))); %节点类型

[angle, l] = cart2pol(node(ele(1:n_ele,2),1)-node(ele(1:n_ele,1),1), node(ele(1:n_ele,2),2)-node(ele(1:n_ele,1),2));
angle = angle*180/pi;

%组装总体刚度矩阵   &   总体质量矩阵
dof=n_node*3+sum(ele(:,3)==2|ele(:,3)==3);             %自由度数，梁单元每个节点有3个自由度，轴向位移、横向位移、面内转角
global di_dof
global se_dof
di_dof=abs(node_type)+3;%节点原始自由度
se_dof=cumsum(di_dof);%累计节点原始自由度
f=ones(dof,1)*1e8;             %外载荷矩阵，整体坐标系下
f_loc=zeros(6,1);         %外载荷矩阵，局部坐标系下
u=ones(dof,1)*1e6;             %位移矩阵
K=zeros(dof);                  %总体刚度矩阵
M=zeros(dof);                   %总体质量矩阵
stress=zeros(n_ele,1);         %单元应力矩阵
node_used=zeros(1,n_node);  %节点使用于铰接或滑动node2的次数
node_used_c=zeros(1,n_node); %copy
for i=1:n_ele
    k_ele=PlaneFrameElementStiffness(A(i),E(i),I(i),l(i),angle(i),ele(i,3));
    K_ele(i)={k_ele};
    [K,node_used(ele(i,2))]=assemPlaneFrame(K,k_ele,ele(i,1),ele(i,2),node_used(ele(i,2)));
    m_ele=PlaneDynamicElementMass(A(i),rou(i),l(i),angle(i));
    M_ele(i)={m_ele};
    [M,node_used_c(ele(i,2))]=assemPlaneFrame(M,m_ele,ele(i,1),ele(i,2),node_used_c(ele(i,2)));
end

%力边界条件（节点力）
f(4)=0000;
f(5)=0;
f(6)=0;
f(7)=0000;
f(8)=0;
f(9)=00000;
f(10)=0;
f(4:end-3)=zeros(dof-6,1);
% f(6)=0;
% f(7)=-10000;
% f(8)=0;
%分布力
f_d=[0 0 1;
    -1600 0 2;
%     0 0 1;
%     0 0 1;
%     -1600 -1600*2/3 2;
%     -1600*2/3 -1600/3 2;
%     -1600/3 0 2;
%     0 0 1;
%     0 0 1;
    0 0 1];
%位移边界条件
u(1)=0;
u(2)=0;
u(3)=0;
u(end-2)=0;
u(end-1)=0;
u(end)=0;

%单元局部坐标系下的载荷转为节点荷载
node_used=zeros(1,n_node);
enf=[];
Dnf=zeros(dof,1);
for i = 1:n_ele
    [enf{i},Dnf,node_used(ele(i,2))] = EquivalentNodeForce(l(i),angle(i),ele(i,1),ele(i,2),ele(i,3), f_d(i,1), f_d(i,2), f_d(i,3) ...
        ,Dnf,node_used(ele(i,2))  ) ;%%%%%%%改
    
end

%求解未知自由度位移
index=[];      %未知自由度的索引
p=[];          %未知自由度对应的节点力矩阵
for i=1:dof
    if u(i)~=0      %划行划列法
        index=[index,i];
        p=[p;f(i)+Dnf(i)];
    end
end
u(index)=K(index,index)\p;    %高斯消去     %第3类element算错            
f=K*u;

%单元局部坐标系下的载荷
node_used(:)=0;
for i=1:n_ele
    %ul=u([3*ele(i,1)-2, 3*ele(i,1)-1, 3*ele(i,1), 3*ele(i,2)-2, 3*ele(i,2)-1, 3*ele(i,2)]);%%%%%%%改
    if  node_type(ele(i,2))==0
        ul=u([se_dof(ele(i,1))-di_dof(ele(i,1))+1:se_dof(ele(i,1))-di_dof(ele(i,1))+3,...
            se_dof(ele(i,2))-di_dof(ele(i,2))+1:se_dof(ele(i,2))-di_dof(ele(i,2))+3]);
    elseif node_type(ele(i,2))~=0
        ul=u([se_dof(ele(i,1))-di_dof(ele(i,1))+1:se_dof(ele(i,1))-di_dof(ele(i,1))+3,...
            se_dof(ele(i,2))-di_dof(ele(i,2))+1:se_dof(ele(i,2))-di_dof(ele(i,2))+2,...
            se_dof(ele(i,2))-di_dof(ele(i,2))+4+node_used(ele(i,2))]);
        node_used(ele(i,2))=node_used(ele(i,2))+1;
    end
    %f_loc=PlaneFrameElementForce(A,E,I(i),l(i),theta(i),ul);   %单元应力局部坐标系下的节点力
    if f_d(i,1)~=0 | f_d(i,2)~=0 
        qEStress(i,:)=stressPlaneFrame(A(i),E(i),I(i),l(i),angle(i),ele(i,3),linspace(0,0,6)',Dnf([3*ele(i,1)-2:3*ele(i,1),3*ele(i,2)-2:3*ele(i,2)]));
        ElementStress(i,:)=stressPlaneFrame(A(i),E(i),I(i),l(i),angle(i),ele(i,3),ul,Dnf([3*ele(i,1)-2:3*ele(i,1),3*ele(i,2)-2:3*ele(i,2)]));
    else
        qEStress(i,:)=linspace(0,0,6);
        ElementStress(i,:)=stressPlaneFrame(A(i),E(i),I(i),l(i),angle(i),ele(i,3),ul,linspace(0,0,6)');
    end
end

%模态
[Xmodal,Omiga2]=eig(K(index,index),M(index,index));
omiga=diag(Omiga2^0.5);
[omiga,sort_omiga]=sort(omiga);
frequency=omiga/2/pi
Xmodal=Xmodal(:,sort_omiga);
XFactor=diag(Xmodal'*M(index,index)*Xmodal);
Xmodaln=Xmodal*inv(sqrt(diag(XFactor))); % Eigenvectors are normalized
%Xmodal=Xmodal(sort_omiga,:);
for i=1:length(index);
    umodal(:,i)=u;
    umodal(index,i)=Xmodal(:,i);
end

%Newmark-beta
C=0.05*M+0.02*K;
f0=100e3;
t1=5;
nt=300;
dt=10^(floor(log10(0.5/max(frequency))));
gama=0.5;
beta=1/6;
a0=1/ beta/dt/dt;
a1=gama/ beta /dt;
a2=1/beta/dt;
a3=1/2/ beta -1;
a4=gama/beta-1;
a5=dt/2*(gama/ beta-2);
a6=dt*(1-gama);
a7=dt*gama;

displacement=zeros(length(index),nt+1);
velocity=zeros(length(index),nt+1);
acceleration=zeros(length(index),nt+1);
tt=zeros(1,nt+1);

for i=2:(nt+1)
    tt(1,1)=1*dt;
    tt(1,i)=i*dt;
    t=(i-1)*dt;
    ft=zeros(length(index),1);
if (t<t1*2)
    ft(1)=f0*sin(4*pi*t/t1);
end
    Ke=K(index,index)+a0*M(index,index)+a1*C(index,index);
    fe=ft+M(index,index)*(a0*displacement(:,i-1)+a2*velocity(:,i-1)+a3*acceleration(:,i-1))+......
        C(index,index)*(a1*displacement(:,i-1)+a4*velocity(:,i-1)+a5*acceleration(:,i-1));
    displacement(:,i)=inv(Ke)*fe;
    acceleration(:,i)=a0*(displacement(:,i)-displacement(:,i-1))-a2*velocity(:,i-1)-a3*acceleration(:,i-1);
    velocity(:,i)=velocity(:,i-1)+a6*acceleration(:,i-1)+a7*acceleration(:,i);
end
udynamic=zeros(dof,nt+1);
udynamic(index,:)=displacement;

%%后处理
nodetj=[]; 
utj=[]; 
paral=zeros(6,n_ele);
gpoint_f=['n','o','d'];
node_used(:)=0;
for i=1:length(ele)
    c=cosd(angle(i));   s=sind(angle(i));
    C1=A(i)*E(i)/l(i);   C2=E(i)*I(i)/l(i)^3;
    T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
    0 0 0 -s c 0;
    0 0 0 0 0 1]; 
    nodetj=T*[0,0,0,node(ele(i,2),1:3)-node(ele(i,1),1:3)]';
    
     %全杆内力
     ftj=T*f([3*ele(i,1)-2, 3*ele(i,1)-1, 3*ele(i,1), 3*ele(i,2)-2, 3*ele(i,2)-1, 3*ele(i,2)]);
     lun0=[nodetj(1).^3, nodetj(1).^2, nodetj(1), 1;
            3*nodetj(1).^2, 2*nodetj(1), 1, 0; 
            nodetj(4).^3, nodetj(4).^2, nodetj(4), 1;
            3*nodetj(4).^2, 2*nodetj(4), 1, 0];
     paraM(:,i)=lun0\[nodetj(2)-ElementStress(i,3),    nodetj(3)-ElementStress(i,2),...
         nodetj(5)+ElementStress(i,6),     nodetj(6)-ElementStress(i,5)]';
    if f_d(i,3)==3
         paraM(:,i)=lun0\[nodetj(2)-ElementStress(i,3),    nodetj(3)-ElementStress(i,2)+f_d(i,1),...
         nodetj(5)+ElementStress(i,6),     nodetj(6)-ElementStress(i,5)+f_d(i,2)]';
     end
     x0tn=nodetj(1):2:nodetj(4);
     lx0tn=[x0tn;x0tn];
     lx0tn=lx0tn(:)';
     M0tn=[paraM(:,i)./5e4]'*[x0tn.^3; x0tn.^2; x0tn; 0*x0tn+1];
     lM0tn=[M0tn;linspace(0,0,length(M0tn))];
     lM0tn=lM0tn(:)';
     lMN=[c,-s;s,c]*[lx0tn;lM0tn]+node(ele(i,1),1:2)';
     lx0n=lMN(1,:); lMn=lMN(2,:);
     
     
     %全杆位移
     node1=ele(i,1);    node2=ele(i,2);    node2_used=node_used(ele(i,2));
     d(1:3)=[se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3];
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
    utj=T*u(d);
    %utj=T*u([3*ele(i,1)-2, 3*ele(i,1)-1, 3*ele(i,1), 3*ele(i,2)-2, 3*ele(i,2)-1, 3*ele(i,2)]);
    xt1=nodetj(1)+utj(1)*20;
    yt1=nodetj(2)+utj(2)*20;
    yxt1=atan(utj(3)*20);
    yxttq1=atan(qEStress(i,3)/E(i)/I(i)*20);
    yxtttq1=atan(qEStress(i,2)/E(i)/I(i)*20);
    xt2=nodetj(4)+utj(4)*20;
    yt2=nodetj(5)+utj(5)*20;
    yxt2=atan(utj(6)*20);
    yxttq2=-atan(qEStress(i,6)/E(i)/I(i)*20);
    yxtttq2=atan(qEStress(i,5)/E(i)/I(i)*20);
    lun=[nodetj(1).^5, nodetj(1).^4, nodetj(1).^3, nodetj(1).^2, nodetj(1), 1; 
                5*nodetj(1).^4, 4*nodetj(1).^3, 3*nodetj(1).^2, 2*nodetj(1), 1, 0;
                20*nodetj(1).^3, 12*nodetj(1).^2, 6*nodetj(1), 2, 0, 0; 
                nodetj(4).^5, nodetj(4).^4, nodetj(4).^3, nodetj(4).^2, nodetj(4), 1;
                5*nodetj(4).^4, 4*nodetj(4).^3, 3*nodetj(4).^2, 2*nodetj(4), 1, 0;
                20*nodetj(4).^3, 12*nodetj(4).^2, 6*nodetj(4), 2, 0, 0];
    switch ele(i,3)
        case 1
        lun=[nodetj(1).^5, nodetj(1).^4, nodetj(1).^3, nodetj(1).^2, nodetj(1), 1; 
                5*nodetj(1).^4, 4*nodetj(1).^3, 3*nodetj(1).^2, 2*nodetj(1), 1, 0;
                20*nodetj(1).^3, 12*nodetj(1).^2, 6*nodetj(1), 2, 0, 0; 
                nodetj(4).^5, nodetj(4).^4, nodetj(4).^3, nodetj(4).^2, nodetj(4), 1;
                5*nodetj(4).^4, 4*nodetj(4).^3, 3*nodetj(4).^2, 2*nodetj(4), 1, 0;
                20*nodetj(4).^3, 12*nodetj(4).^2, 6*nodetj(4), 2, 0, 0];
        
         paral(3:6,i)=lun0\[yt1,yxt1,yt2,yxt2]';
         paralq(:,i)=lun\[0,0,yxttq1,0,0,yxttq2]';
         paral(:,i)=paral(:,i)+paralq(:,i);
%         paral(:,i)=[paraM(:,i).*[1/20/E(i)/I(i),1/12/E(i)/I(i),1/6/E(i)/I(i),1/2/E(i)/I(i)]';atan(utj(2));utj(1)];
         
        case 2
        lun=[nodetj(1).^5, nodetj(1).^4, nodetj(1).^3, nodetj(1).^2, nodetj(1), 1; 
                5*nodetj(1).^4, 4*nodetj(1).^3, 3*nodetj(1).^2, 2*nodetj(1), 1, 0;
                20*nodetj(1).^3, 12*nodetj(1).^2, 6*nodetj(1), 2, 0, 0; 
                nodetj(4).^5, nodetj(4).^4, nodetj(4).^3, nodetj(4).^2, nodetj(4), 1;
                20*nodetj(4).^3, 12*nodetj(4).^2, 6*nodetj(4), 2, 0, 0;
                60*nodetj(4).^2, 24*nodetj(4), 6, 0, 0, 0];
         paral(4:6,i)=lun0(1:3,2:4)\[yt1,yxt1,yt2]';
         paralq(:,i)=lun\[0,0,yxttq1,0,yxttq2,yxtttq2]';
         paral(:,i)=paral(:,i)+paralq(:,i);

%          paral(:,i)=[paraM(:,i).*[1/20/E(i)/I(i),1/12/E(i)/I(i),1/6/E(i)/I(i),1/2/E(i)/I(i)]';atan(utj(2));utj(1)];
         
        case 3
        lun=[nodetj(1).^5, nodetj(1).^4, nodetj(1).^3, nodetj(1).^2, nodetj(1), 1; 
                5*nodetj(1).^4, 4*nodetj(1).^3, 3*nodetj(1).^2, 2*nodetj(1), 1, 0;
                20*nodetj(1).^3, 12*nodetj(1).^2, 6*nodetj(1), 2, 0, 0; 
                5*nodetj(4).^4, 4*nodetj(4).^3, 3*nodetj(4).^2, 2*nodetj(4), 1, 0;
                20*nodetj(4).^3, 12*nodetj(4).^2, 6*nodetj(4), 2, 0, 0;
                60*nodetj(4).^2, 24*nodetj(4), 6, 0, 0, 0];
             paral(4:6,i)=lun0([1:2,4],2:4)\[yt1,yxt1,yxt2]';
             paralq(:,i)=lun\[0,0,yxttq1,0,0,yxttq2]';
             paral(:,i)=paral(:,i)+paralq(:,i);
%          paral(:,i)=[paraM(:,i).*[1/20/E(i)/I(i),1/12/E(i)/I(i),1/6/E(i)/I(i),1/2/E(i)/I(i)]';atan(utj(2));utj(1)];
         
        case 0
%         lun=[xt1, 1; 
%                 xt2, 1];
%          paral(i)={lun\[yt1,yt2]'};
%          xtn=xt1:1:xt2;
%          ytn=[cell2mat(paral(i))]'*[xtn; 0*xtn+1];
            lun=[nodetj(1).^5, nodetj(1).^4, nodetj(1).^3, nodetj(1).^2, nodetj(1), 1; 
                5*nodetj(1).^4, 4*nodetj(1).^3, 3*nodetj(1).^2, 2*nodetj(1), 1, 0;
                60*nodetj(1).^2, 24*nodetj(1), 6, 0, 0, 0; 
                nodetj(4).^5, nodetj(4).^4, nodetj(4).^3, nodetj(4).^2, nodetj(4), 1;
                5*nodetj(4).^4, 4*nodetj(4).^3, 3*nodetj(4).^2, 2*nodetj(4), 1, 0;
                60*nodetj(4).^2, 24*nodetj(4), 6, 0, 0, 0];
            paral(5:6,i)=lun0([1,3],3:4)\[yt1,yt2]';
             paralq(:,i)=lun\[0,0,yxtttq1,0,0,yxtttq2]';
             paral(:,i)=paral(:,i)+paralq(:,i);
    end
    xtn=linspace(xt1,xt2,length(x0tn));
    ytn=[paral(:,i)]'*[x0tn.^5; x0tn.^4; x0tn.^3; x0tn.^2; x0tn; 0*x0tn+1];
    
     XYN=[c,-s;s,c]*[xtn;ytn]+node(ele(i,1),1:2)';
     xn=XYN(1,:); yn=XYN(2,:);
     
     figure(1)%位移图
     %     line(node(ele(i,1:2),1)+u(3*ele(i,1:2)-2)*100,node(ele(i,1:2),2)+u(3*ele(i,1:2)-1)*100,...
     %         'LineWidth',2,'color','b','LineStyle','--','Marker','*')
     if ele(i,3)~=0
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r')
         line(node(ele(i,2),1),node(ele(i,2),2),'LineWidth',2,'color','r','Marker',gpoint_f(ele(i,3)))
     else
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r','Marker',gpoint_f(2))
     end
     
     line(xn,yn, 'LineWidth',2,'color','b','LineStyle','--')
     
     figure(2)%弯矩图
     if ele(i,3)~=0
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r')
         line(node(ele(i,2),1),node(ele(i,2),2),'LineWidth',2,'color','r','Marker',gpoint_f(ele(i,3)))
     else
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r','Marker',gpoint_f(2))
     end
     
     for ii = 1:length(x0tn)
         line(lx0n(2*ii-1:2*ii),lMn(2*ii-1:2*ii),'LineWidth',2,'LineStyle','-')
     end
     
     figure(3)%N阶模态
     N_mo=1;
     [Nxdl1,Nxdl2]=ShapeStrainFuncMatrix(0:0.02:1,l(i));
     utj_mo1=T*umodal(d,N_mo).*1e4;
     xtn_mo1=[Nxdl1*utj_mo1]'+[0:0.02:1].*l(i);
     ytn_mo1=[Nxdl2*utj_mo1]';
     XYN_mo1=[c,-s;s,c]*[xtn_mo1;ytn_mo1]+node(ele(i,1),1:2)';
     xn_mo1=XYN_mo1(1,:); yn_mo1=XYN_mo1(2,:);
     if ele(i,3)~=0
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r')
         line(node(ele(i,2),1),node(ele(i,2),2),'LineWidth',2,'color','r','Marker',gpoint_f(ele(i,3)))
     else
         line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r','Marker',gpoint_f(2))
     end
     line(xn_mo1,yn_mo1, 'LineWidth',2,'color','b','LineStyle','--')

end
figure(1)
title('位移图*20')
axis equal
grid on
figure(2)
title('弯矩图/5e4')
axis equal
grid on
figure(3)
title(['第',mat2str(N_mo),'阶模态'])
text(-30,10,['frequency=',mat2str(frequency(N_mo))])
axis equal
grid on
figure(4)%Newmark-beta
nri_dof=1;
plot(tt(nri_dof,:),displacement(nri_dof,:))
ylabel('位移/[m]','fontsize',12);
grid on;
xlabel('时间/[s]','fontsize',12);
title(['第',mat2str(nri_dof),'实际自由度位移响应时程'],'fontsize',12);

figure(5)
for ntt=0:nt
    
    for i=1:length(ele)
        c=cosd(angle(i));   s=sind(angle(i));
        C1=A(i)*E(i)/l(i);   C2=E(i)*I(i)/l(i)^3;
        T=[ c s 0 0 0 0;    -s c 0 0 0 0;   0 0 1 0 0 0;
            0 0 0 c s 0;    0 0 0 -s c 0;     0 0 0 0 0 1]; 
        node1=ele(i,1);    node2=ele(i,2);    node2_used=node_used(ele(i,2));
         d(1:3)=[se_dof(node1)-di_dof(node1)+1:se_dof(node1)-di_dof(node1)+3];
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
        nodetj=T*[0,0,0,node(ele(i,2),1:3)-node(ele(i,1),1:3)]';
        utj_dynt=T*udynamic(d,ntt+1);
        utj_dynt_mul=utj_dynt.*20;
        xtn_dynt=[Nxdl1*utj_dynt_mul]'+[0:0.02:1].*l(i);
        ytn_dynt=[Nxdl2*utj_dynt_mul]';
        XYN_dynt=[c,-s;s,c]*[xtn_dynt;ytn_dynt]+node(ele(i,1),1:2)';
        xn_dynt=XYN_dynt(1,:); yn_dynt=XYN_dynt(2,:);
        if ele(i,3)~=0
            line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r')
            line(node(ele(i,2),1),node(ele(i,2),2),'LineWidth',2,'color','r','Marker',gpoint_f(ele(i,3)))
        else
            line(node(ele(i,1:2),1),node(ele(i,1:2),2),'LineWidth',2,'color','r','Marker',gpoint_f(2))
        end
        line(xn_dynt,yn_dynt, 'LineWidth',2,'color','b','LineStyle','--')
    end
    axis equal
    grid on
    title(['位移响应*20,t=',mat2str(ntt*dt),'s'])
    pause(dt/10);
    clf(figure(5))
end

ElementStress