
%A simple FEM implementation for the solution to the heat equation in
%a triangular domain. Neumann and Dirichlet BC's are used.

clear all;
close all;
 
h=0.288675;
%N_e: number of elements
%N_p: number of nodes, including boundary nodes
%N_b: number of boundary nodes of Dirichlet type
N_e=16;
N_p=15;
N_b=12;
Nodes=[-0.57735 0;
   -0.288675 0;
   0 0;
   0.288675 0;
   0.57735 0;
   
   -0.433013 0.25;
   -0.144338 0.25;
   0.144338 0.25;
   0.433013 0.25;
   
   -0.288675 0.5;
   0 0.5;
   0.288675 0.5;
   
   -0.144338 0.75;
   0.144338 0.75;
   
   0 1
   ];
%B is the connectivity matrix
B=[1 2 6;
    2 7 6;
    2 3 7;
    3 8 7;
    3 4 8;
    4 9 8;
    4 5 9;
    6 7 10;
    7 11 10;
    7 8 11;
    8 12 11;
    8 9 12;
    10 11 13;
    11 14 13;
    11 12 14;
    13 14 15
    ];
 
Dirichlet_nodes=[1;
                 2;
                 3;
                 4;
                 5;
                 6;
                 9;
                 10;
                 12;
                 13;
                 14;
                 15
                 
                 
                 
                 ];
A_element=[1/sqrt(3) -1/(2*sqrt(3)) -1/(2*sqrt(3));
      -1/(2*sqrt(3)) 1/sqrt(3) -1/(2*sqrt(3));
      -1/(2*sqrt(3)) -1/(2*sqrt(3)) 1/sqrt(3)];
b_element=[h^2/(4*sqrt(3)) h^2/(4*sqrt(3)) h^2/(4*sqrt(3))];
 
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
        i_dirichlet=i_dirichlet+1;
    else
        map(i_map)=i;
        i_map=i_map+1;
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
c_sub=A_sub\b_sub;
c_global=zeros(N_p,1);
for i=1:N_p-N_b
    i0=map(i);
    c_global(i0,1)=c_sub(i,1);
end
%----------------------------------------------------------------------------
 
%postprocessing
%plot the value at nodes 1 to 5
%plot(Nodes(1:5,1),c_global(1:5,1),'-o');
xlabel('x');
ylabel('u');
legend('FEM');
 
%plot the distribution on each element
patch('Vertices',Nodes,'Faces',B,'FaceVertexCData',c_global,'FaceColor','interp');
 
c_global(7)
c_global(8)
c_global(11)

