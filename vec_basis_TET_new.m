clc;clear all;
%% paramters of tetrahedron
La=4;
pos(1,1:3)=[-La/2,-tan(1/6*pi)*La/2,0];
pos(2,1:3)=[ La/2,-tan(1/6*pi)*La/2,0];
pos(3,1:3)=[0,(tan(1/3*pi)-tan(1/6*pi))*La/2,   0];
pos(4,1:3)=[0,0,sqrt(6)*La/3];
F = [1, 3, 2 ;  4, 2 ,3;  1 ,4, 3;1, 2 ,4] ;

for m=1:4
x(m)=pos(m,1);
y(m)=pos(m,2);
z(m)=pos(m,3);
end
aa(1)=x(2)*(y(3)*z(4)-z(3)*y(4))-x(3)*(y(2)*z(4)-z(2)*y(4))+x(4)*(y(2)*z(3)-z(2)*y(3));
aa(2)=-x(1)*(y(3)*z(4)-z(3)*y(4))+x(3)*(y(1)*z(4)-z(1)*y(4))-x(4)*(y(1)*z(3)-z(1)*y(3));
aa(3)=x(1)*(y(2)*z(4)-z(2)*y(4))-x(2)*(y(1)*z(4)-z(1)*y(4))+x(4)*(y(1)*z(2)-z(1)*y(2));
aa(4)=-x(1)*(y(2)*z(3)-z(2)*y(3))+x(2)*(y(1)*z(3)-z(1)*y(3))-x(3)*(y(1)*z(2)-z(1)*y(2));

deltae=(aa(1)+aa(2)+aa(3)+aa(4))/6;

bb(1)=-(y(3)*z(4)-z(3)*y(4))+(y(2)*z(4)-z(2)*y(4))-(y(2)*z(3)-z(2)*y(3));
bb(2)=(y(3)*z(4)-z(3)*y(4))-(y(1)*z(4)-z(1)*y(4))+(y(1)*z(3)-z(1)*y(3));
bb(3)=-(y(2)*z(4)-z(2)*y(4))+(y(1)*z(4)-z(1)*y(4))-(y(1)*z(2)-z(1)*y(2));
bb(4)=(y(2)*z(3)-z(2)*y(3))-(y(1)*z(3)-z(1)*y(3))+(y(1)*z(2)-z(1)*y(2));

cc(1)=(x(3)*z(4)-z(3)*x(4))-(x(2)*z(4)-z(2)*x(4))+(x(2)*z(3)-z(2)*x(3));
cc(2)=-(x(3)*z(4)-z(3)*x(4))+(x(1)*z(4)-z(1)*x(4))-(x(1)*z(3)-z(1)*x(3));
cc(3)=(x(2)*z(4)-z(2)*x(4))-(x(1)*z(4)-z(1)*x(4))+(x(1)*z(2)-z(1)*x(2));
cc(4)=-(x(2)*z(3)-z(2)*x(3))+(x(1)*z(3)-z(1)*x(3))-(x(1)*z(2)-z(1)*x(2));

dd(1)=-(x(3)*y(4)-y(3)*x(4))+(x(2)*y(4)-y(2)*x(4))-(x(2)*y(3)-y(2)*x(3));
dd(2)=(x(3)*y(4)-y(3)*x(4))-(x(1)*y(4)-y(1)*x(4))+(x(1)*y(3)-y(1)*x(3));
dd(3)=-(x(2)*y(4)-y(2)*x(4))+(x(1)*y(4)-y(1)*x(4))-(x(1)*y(2)-y(1)*x(2));
dd(4)=(x(2)*y(3)-y(2)*x(3))-(x(1)*y(3)-y(1)*x(3))+(x(1)*y(2)-y(1)*x(2));
%% Sampling point in this tetrahedron
TRI=1;
N_P=5;nn=0;
i1=1;i2=3;      %It is mass if the basis functions are all drawed, so I choose one of them.
le=sqrt((pos(i1,1)-pos(i2,1))^2+(pos(i1,2)-pos(i2,2))^2+(pos(i1,3)-pos(i2,3))^2);   %length of this basis function
for ss=1:4
j1=F(ss,1);j2=F(ss,2);j3=F(ss,3);
    for m=0:N_P+1
        for n=0:N_P+1
            if((1-n/(N_P+1)-m/(N_P+1))>=-0.01)
            nn=nn+1;
            pos_t(nn,1)=pos(j1,1)*(m/(N_P+1))+pos(j2,1)*(n/(N_P+1))+pos(j3,1)*(1-n/(N_P+1)-m/(N_P+1));
            pos_t(nn,2)=pos(j1,2)*(m/(N_P+1))+pos(j2,2)*(n/(N_P+1))+pos(j3,2)*(1-n/(N_P+1)-m/(N_P+1));
            pos_t(nn,3)=pos(j1,3)*(m/(N_P+1))+pos(j2,3)*(n/(N_P+1))+pos(j3,3)*(1-n/(N_P+1)-m/(N_P+1));
            end
            x_c=pos_t(nn,1);y_c=pos_t(nn,2);z_c=pos_t(nn,3);
            tmp1=aa(i1)+bb(i1)*x_c+cc(i1)*y_c+dd(i1)*z_c;
            tmp2=aa(i2)+bb(i2)*x_c+cc(i2)*y_c+dd(i2)*z_c;
            Nx(nn)=le/(6*deltae)^2*(tmp1*bb(i2)-tmp2*bb(i1));
            Ny(nn)=le/(6*deltae)^2*(tmp1*cc(i2)-tmp2*cc(i1));
            Nz(nn)=le/(6*deltae)^2*(tmp1*dd(i2)-tmp2*dd(i1));
        end
    end
end
nn=nn+1;
pos_t(nn,1)=sum(pos(:,1))/4;
pos_t(nn,2)=sum(pos(:,2))/4;
pos_t(nn,3)=sum(pos(:,3))/4;
x_c=pos_t(nn,1);y_c=pos_t(nn,2);z_c=pos_t(nn,3);
tmp1=aa(i1)+bb(i1)*x_c+cc(i1)*y_c+dd(i1)*z_c;
tmp2=aa(i2)+bb(i2)*x_c+cc(i2)*y_c+dd(i2)*z_c;
Nx(nn)=le/(6*deltae)^2*(tmp1*bb(i2)-tmp2*bb(i1));
Ny(nn)=le/(6*deltae)^2*(tmp1*cc(i2)-tmp2*cc(i1));
Nz(nn)=le/(6*deltae)^2*(tmp1*dd(i2)-tmp2*dd(i1));
%% drawe
h1=figure(1);
aaa=patch('Faces',F,'Vertices',pos,  'EdgeColor','R','facealpha',0.4);%,'LineStyle','none');
xlabel('{\it\bfX}{\rm\bf/m}','fontsize',28,'fontweight','b','FontName','Times New Roman');
ylabel('{\it\bfY}{\rm\bf/m}','fontsize',28,'fontweight','b','FontName','Times New Roman');
zlabel('{\it\bfZ}{\rm\bf/m}','fontsize',28,'fontweight','b','FontName','Times New Roman');
set(gcf, 'renderer', 'zbuffer');
set(gcf,'color','white')
axis off;grid on;
axis equal;
alpha(0.2)
hold on
scatter3(pos_t(:,1),pos_t(:,2),pos_t(:,3),'k');
L=1:size(pos,1);
scatter3(pos(:,1),pos(:,2),pos(:,3),210,'filled');                   
text(pos(:,1),pos(:,2),pos(:,3),num2str(L(:)),'fontsize',28)         
hcone = coneplot(pos_t(:,1),pos_t(:,2),pos_t(:,3),Nx(:),Ny(:),Nz(:),0.1,'nointerp') ;
hcone.FaceColor = [.2 .7 1];
hcone.EdgeColor = 'none';
plot3 (pos([i1 i2],1),pos([i1 i2],2),pos([i1 i2],3),'color','G','linewidth',6)       
camlight, lighting gouraud
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	 
jFrame = get(h1,'JavaFrame');	 
pause(0.1);					 
set(jFrame,'Maximized',1);	 
pause(0.1);					 
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		 
text(pos_t(nn,1),pos_t(nn,2),pos_t(nn,3),'Q. Yang (Xidian)','fontsize',24,'fontweight','b','FontName','Times New Roman')        
for i=1:2:360
if ((i)~=0)&(mod((i),10)==0)
fprintf('%5i, %5i\n',i,len)
end
view(i,10)      %output GIF
GIF(i)=getframe(gcf);
    im = frame2im(GIF(i));
    [I,map] = rgb2ind(im,256);
    %GIF 
    if i == 1
        imwrite(I,map,'test.gif','GIF', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map,'test.gif','GIF','WriteMode','append','DelayTime',0.1);
    end    
end