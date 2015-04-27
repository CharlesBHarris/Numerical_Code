%Math 822 Final Project, Charles Harris, Finite Element Method, A=4




%h=1/4;
%N_e: number of elements
%N_p: number of nodes, including boundary nodes
%N_b: number of boundary nodes of Dirichlet type
N_e=96;
N_p=65;
N_b=5;

A=4;

Nodes=zeros(N_p,2);
for i=1:A
 for j=1:A+1
    Nodes(5*(i-1)+j,1)= (j-1)*(2/A);
    Nodes(5*(i-1)+j,2)= (i-1)*(2/A)-2;
 end
end

for i=1:A+1
 for j=1:(2*A+1)
    Nodes(20+9*(i-1)+j,1)= (j-1)*(2/A)-2;
    Nodes(20+9*(i-1)+j,2)= (i-1)*(2/A);
 end
end


%B is the connectivity matrix
B=[1 6 2;
    2 6 7;
    2 7 3;
    3 7 8;
    3 8 4;
    4 8 9;
    4 9 5;
    5 9 10;
    6 11 7;
    7 11 12;
    7 12 8;
    8 12 13;
    8 13 9;
    9 13 14;
    9 14 10;
    10 14 15;
    11 16 12;
    12 16 17;
    12 17 13;
    13 17 18;
    13 18 14;
    14 18 19;
    14 19 15;
    15 19 20;
    16 25 17;
    17 25 26;
    17 26 18;
    18 26 27;
    18 27 19;
    19 27 28;
    19 28 20;
    20 28 29;
    21 30 22;
    22 30 31;
    22 31 23;
    23 31 32;
    23 32 24;
    24 32 33;
    24 33 25;
    25 33 34;
    25 34 26;
    26 34 35;
    26 35 27;
    27 35 36;
    27 36 28;
    28 36 37;
    28 37 29;
    29 37 38;
    30 39 31;
    31 39 40;
    31 40 32;
    32 40 41;
    32 41 33;
    33 41 42;
    33 42 34;
    34 42 43;
    34 43 35;
    35 43 44;
    35 44 36;
    36 44 45;
    36 45 37;
    37 45 46;
    37 46 38;
    38 46 47;
    39 48 40;
    40 48 49;
    40 49 41;
    41 49 50;
    41 50 42;
    42 50 51;
    42 51 43;
    43 51 52;
    43 52 44;
    44 52 53;
    44 53 45;
    45 53 54;
    45 54 46;
    46 54 55;
    46 55 47;
    47 55 56;
    48 57 49;
    49 57 58;
    49 58 50;
    50 58 59;
    50 59 51;
    51 59 60;
    51 60 52;
    52 60 61;
    52 61 53;
    53 61 62;
    53 62 54;
    54 62 63;
    54 63 55;
    55 63 64;
    55 64 56;
    56 64 65];

Dirichlet_nodes=[21;
                 30;
                 39;
                 48;
                 57];
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
for i=1:5;
    j=Dirichlet_nodes(i);
    C(:,i)=A_global(:,j);
end
C(21,:)=[];
C(29,:)=[];
C(37,:)=[];
C(45,:)=[];
C(53,:)=[];
E=C(:,2)*(2*.5-.5^2)+C(:,3)+C(:,4)*(2*1.5-1.5^2);
e=-E;


c_sub=A_sub\e;
c_new=[c_sub(1:20);0;c_sub(21:28);0.75;c_sub(29:36);1;c_sub(37:44);0.75;c_sub(45:52);0;c_sub(53:60)];
disp('soln')
c_new(35)
c_new(45)
c_new(55)





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





