%This code implements a simple FEM method for calculating the Navier-Stokes
%flow of a hypothetical fluid in an irregular rectangular domain. A Neumann
%boundary conditions is utilized.



%h=1/4;
%N_e: number of elements
%N_p: number of nodes, including boundary nodes
%N_b: number of boundary nodes of Dirichlet type
N_e=384;
N_p=225;
N_b=9;

A=8;

Nodes=zeros(N_p,2);
for i=1:A
 for j=1:A+1
    Nodes(9*(i-1)+j,1)= (j-1)*(2/A);
    Nodes(9*(i-1)+j,2)= (i-1)*(2/A)-2;
 end
end

for i=1:A+1
 for j=1:(2*A+1)
    Nodes(72+17*(i-1)+j,1)= (j-1)*(2/A)-2;
    Nodes(72+17*(i-1)+j,2)= (i-1)*(2/A);
 end
end



W_1=zeros(16,3);
for i=1:16
   if mod(i,2)==0
       W_1(i,1)=1+(i/2);
       W_1(i,2)=10+(i/2);
       W_1(i,3)=9+(i/2);
   else
       W_1(i,1)=((i+1)/2);
       W_1(i,2)=1+((1+i)/2);
       W_1(i,3)=9+((1+i)/2);
   end
end


W_2=zeros(32,3);
for i=1:32
   if mod(i,2)==0
       W_2(i,1)=73+(i/2);
       W_2(i,2)=90+(i/2);
       W_2(i,3)=89+(i/2);
   else
       W_2(i,1)=72+((i+1)/2);
       W_2(i,2)=73+((1+i)/2);
       W_2(i,3)=89+((1+i)/2);
   end
end

W_3=zeros(16,3);
for i=1:16
   if mod(i,2)==0
       W_3(i,1)=64+(i/2);
       W_3(i,2)=81+(i/2);
       W_3(i,3)=80+(i/2);
   else
       W_3(i,1)=63+((i+1)/2);
       W_3(i,2)=64+((1+i)/2);
       W_3(i,3)=80+((1+i)/2);
   end
end



B=zeros(384,3);
for i=1:7
    for k=1:16
    for j=1:3     
   B(k+16*(i-1),j)=W_1(k,j)+9*(i-1);
    end
    end
end

for k=1:16
    for j=1:3
    B(112+k,j)=W_3(k,j);  
        
    end
end

for i=1:8
    for k=1:32
    for j=1:3     
   B(k+32*(i-1)+128,j)=W_2(k,j)+17*(i-1);
    end
    end
end

B






%];

Dirichlet_nodes=[73;
                 90;
                 107;
                 124;
                 141;
                 158;
                 175;
                 192;
                 209];
A_element=[1 -1/2 -1/2;
      -1/2 1/2 0;
      -1/2 0 1/2];
b_element=[0 0 0];

%----assemble the global matrix and solution is stored in c_global----------------
A_global=zeros(N_p);
b_global=zeros(N_p,1);

for i_e=1:N_e
    for m=1:3
        for n=1:3
            i=B(i_e,m);
            j=B(i_e,n);
            A_global(i,j)=A_global(i,j)+A_element(m,n);
        end
    end
    for m=1:3
        i=B(i_e,m);
        b_global(i)=b_global(i)+b_element(m);
    end
end

map=(1:(N_p-N_b))';
i_map=1;
i_dirichlet=1;
for i=1:N_p
    if(i==Dirichlet_nodes(i_dirichlet))
        if i_dirichlet < N_b
        i_dirichlet=i_dirichlet+1
        end
    else
        map(i_map)=i;
        i_map=i_map+1
    end
end
A_sub=zeros(N_p-N_b);
b_sub=zeros(N_p-N_b,1);
for i=1:N_p-N_b
    for j=1:N_p-N_b
        i0=map(i);
        j0=map(j);
        A_sub(i,j)=A_global(i0,j0);
    end
end
for i=1:N_p-N_b
    i0=map(i);
    b_sub(i,1)=b_global(i0,1);
end
%c_sub=A_sub\b_sub;

C=zeros(N_p,5);
for i=1:9;
    j=Dirichlet_nodes(i);
    C(:,i)=A_global(:,j);
end
C(73,:)=[];
C(89,:)=[];
C(105,:)=[];
C(121,:)=[];
C(137,:)=[];
C(153,:)=[];
C(169,:)=[];
C(185,:)=[];
C(201,:)=[];

E=C(:,2)*(2*.25-.25^2)+C(:,3)*(2*.5-.5^2)+C(:,4)*(2*.75-.75^2)+C(:,5)+C(:,6)*(2*1.25-1.25^2)+C(:,7)*(2*1.5-1.5^2)+C(:,8)*(2*1.75-1.75^2);
e=-E;


c_sub=A_sub\e;
c_new=[c_sub(1:72);0;c_sub(73:88);0.4375;c_sub(89:104);.75;c_sub(105:120);0.9375;c_sub(121:136);1;c_sub(137:152);.9375;c_sub(153:168);.75;c_sub(169:184);0.4375;c_sub(185:200);0;c_sub(201:216)];
disp('soln')
c_new(117)
c_new(153)
c_new(189)




% c_global=zeros(N_p,1);
% for i=1:N_p-N_b
%     i0=map(i);
%     c_global(i0,1)=c_sub(i,1);
% end
% %----------------------------------------------------------------------------
% 
% %postprocessing
% %plot the value at nodes 1 to 5
% exact=ones(5,1)*0.5;
% for n=1:20
%     alpha=(2*n-1)*pi/2;
%     exact=exact+2*(-1)^n*cosh(alpha*Nodes(1:5,1))/(alpha^3)/cosh(alpha);
% end
% exact
% plot(Nodes(1:5,1),c_global(1:5,1),'-o',Nodes(1:5,1),exact,'-s');
% xlabel('x');
% ylabel('u');
legend('FEM');

%plot the distribution on each element
patch('Vertices',Nodes,'Faces',B,'FaceVertexCData',c_new,'FaceColor','interp');





