%This code implements a simple FEM method for calculating the Navier-Stokes
%flow of a hypothetical fluid in an irregular rectangdomain. A mixture of
%Dirchelet and Neumann boundary conditions are utilized.



%------------INPUT GEOMETRY AND BOUNDRY CONDITIONS
clear;
N_x=8;
N_nodes=8*N_x;
h=2/(N_x);

%coordinates of element endpoints
x_elements=[];
for j=1:N_x
    x_elements=[x_elements; (j-1)*h -2];
end
for j=1:(2*N_x+1)
    x_elements=[x_elements; 2 (-N_x+j-1)*h];
end
for j=(2*N_x):-1:1
    x_elements=[x_elements; (-N_x+j-1)*h 2];
end
for j=(N_x-1):-1:1
    x_elements=[x_elements; -2 j*h];
end
for j=N_x:-1:1
    x_elements=[x_elements; -j*h 0];
end
for j=1:(N_x)
    x_elements=[x_elements; 0 (-j+1)*h];
end

x_elements=[x_elements;0 -2];

%specify boundary conditions
Dirichlet_BC=[];
for j=1:9
    z=(N_x-(j-1))*h;
    Dirichlet_BC=[Dirichlet_BC; 40+j 2*z-z^2];
end

N_Dirichlet=N_x+1;
Newmann_BC=[];
for j=1:55
   
    if j>40;
        
    j=j+9;
    
    end
    Newmann_BC=[Newmann_BC; j 0];
end
N_Newmann=N_nodes-(N_x+1);

%------SOLVE THE BOUNDARY INTEGRAL EQUATION-------------------------

%coordinates of nodal points
x_nodes=[];
for j=1:N_nodes
    x_nodes=[x_nodes; 0.5*(x_elements(j,1:2)+x_elements(j+1,1:2))];
end

%form normal vectors
n_vector=[];
L=(1:N_nodes)*0;
for j=1:N_nodes
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    L(j)=norm(beta);
    n_vector=[n_vector; beta(2)/L(j) -beta(1)/L(j)];
end


u=zeros(N_nodes,1);
q=zeros(N_nodes,1);

for j=1:N_Dirichlet
    u(Dirichlet_BC(j,1),1)=Dirichlet_BC(j,2);
end
for j=1:N_Newmann
    q(Newmann_BC(j,1),1)=Newmann_BC(j,2);
end

%choose a quadrature
N_p=4;
w=[0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];
p=[-0.861136311594053 -0.339981043584856 0.339981043584856 0.861136311594053];
w_new=w/2;
p_new=(p+1)/2;
tmp2=0*(1:N_p);
%form G and H matrices
for j=1:N_nodes
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    for i=1:N_nodes
        if(i ~= j)
            alpha=x_elements(j,1:2)-x_nodes(i,1:2);
            
            %G matrix
            for k=1:N_p
                tmp2(k)=log(norm(alpha+beta*p_new(k)));
            end
            G(i,j)=-(L(j)/2/pi)*(w_new*tmp2');
            
            %H matrix
            for k=1:N_p
                tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
            end
            H(i,j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
            
        else
            G(i,i)=-(L(j)/2/pi)*(-1+log(L(j)/2));
            H(i,i)=0.5;
        end
    end
end


for j=1:N_Dirichlet
    jd=Dirichlet_BC(j,1);
    q(jd,1)=u(jd,1);
    for i=1:N_nodes
        swap=H(i,jd);
        H(i,jd)=-G(i,jd);
        G(i,jd)=-swap;
    end
end
b=G*q;
u=H\b;
for j=1:N_Dirichlet
    jd=Dirichlet_BC(j,1);
    swap=u(jd,1);
    u(jd,1)=q(jd,1);
    q(jd,1)=swap;
end

%---------POST PROCESSING-----------------------------------------------

%choose the points for output of the solution
nx=3;
x_output=[.5 1 1.5];
y_output=[.5 1 1.5];
u_output=0*(1:nx);


for i=1:(nx)
    xp=x_output(i);
    yp=y_output(i);
    for j=1:N_nodes
        beta=x_elements(j+1,1:2)-x_elements(j,1:2);
        alpha=x_elements(j,1:2)-[xp yp];
        %g vector
        for k=1:N_p
            tmp2(k)=log(norm(alpha+beta*p_new(k)));
        end
        g_vector(j)=-(L(j)/2/pi)*(w_new*tmp2');
            
        %H matrix
        for k=1:N_p
            tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
        end
        h_vector(j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
    end
    u_output(i)=g_vector*q-h_vector*u;
end
u_output
% norm(u_output-x_output)
% plot(x_output,u_output,'+',x_output,x_output,'--');
% legend('numerical','exact',4);
%plot(x_output,u_output-x_output);

%for contour plots
vortices=[];

for i=1:N_nodes
    vortices=[vortices; x_elements(i,1:2)];
end
for i=1:161;
   if i<=7;
       vortices=[vortices; i*h -1.75];
   elseif (i>7) & (i<=14)
    vortices=[vortices; (i-7)*h -1.5];
    elseif (i>14) & (i<=21)
    vortices=[vortices; (i-14)*h -1.25];
    elseif (i>21) & (i<=28)
    vortices=[vortices; (i-21)*h -1];
    elseif (i>28) & (i<=35)
    vortices=[vortices; (i-28)*h -0.75];
    elseif (i>35) & (i<=42)
    vortices=[vortices; (i-35)*h -0.5];
    elseif (i>42) & (i<=49)
    vortices=[vortices; (i-42)*h -0.25];
    elseif (i>49) & (i<=56)
    vortices=[vortices; (i-49)*h 0];
    elseif (i>56) & (i<=71)
    vortices=[vortices; (-8+(i-56))*h 0.25];
    elseif (i>71) & (i<=86)
    vortices=[vortices; (-8+(i-71))*h 0.5];
    elseif (i>86) & (i<=101)
    vortices=[vortices; (-8+(i-86))*h 0.75];
    elseif (i>101) & (i<=116)
    vortices=[vortices; (-8+(i-101))*h 1];
    elseif (i>116) & (i<=131)
    vortices=[vortices; (-8+(i-116))*h 1.25];
    elseif (i>131) & (i<=146)
    vortices=[vortices; (-8+(i-131))*h 1.5];
    elseif (i>146) & (i<=161)
    vortices=[vortices; (-8+(i-146))*h 1.75];
   end
 end









B=[%1
    1 2 64;
    2 65 64;
    2 3 65;
    3 66 65;
    3 4 66;
    4 67 66;
    4 5 67;
    5 68 67;
    5 6 68;
    6 69 68;
    6 7 69;
    7 70 69;
    7 8 70;
    8 71 70;
    8 9 71;
    9 10 71;
    %2
    64 65 63;
    65 72 63;
    65 66 72;
    66 73 72;
    66 67 73;
    67 74 73;
    67 68 74;
    68 75 74;
    68 69 75;
    69 76 75;
    69 70 76;
    70 77 76;
    70 71 77;
    71 78 77;
    71 10 78;
    10 11 78;
    %3
    63 72 62;
    72 79 62;
    72 73 79;
    73 80 79;
    73 74 80;
    74 81 80;
    74 75 81;
    75 82 81;
    75 76 82;
    76 83 82;
    76 77 83;
    77 84 83;
    77 78 84;
    78 85 84;
    78 11 85;
    11 12 85;
    %4
    62 79 61;
    79 86 61;
    79 80 86;
    80 87 86;
    80 81 87;
    81 88 87;
    81 82 88;
    82 89 88;
    82 83 89;
    83 90 89;
    83 84 90;
    84 91 90;
    84 85 91;
    85 92 91;
    85 12 92;
    12 13 92;
    %5
    61 86 60;
    86 93 60;
    86 87 93;
    87 94 93;
    87 88 94;
    88 95 94;
    88 89 95;
    89 96 95;
    89 90 96;
    90 97 96;
    90 91 97;
    91 98 97;
    91 92 98;
    92 99 98;
    92 13 99;
    13 14 99;
    %6
    60 93 59;
    93 100 59;
    93 94 100;
    94 101 100;
    94 95 101;
    95 102 101;
    95 96 102;
    96 103 102;
    96 97 103;
    97 104 103;
    97 98 104;
    98 105 104;
    98 99 105;
    99 106 105;
    99 14 106;
    14 15 106;
    %7
    59 100 58;
    100 107 58;
    100 101 107;
    101 108 107;
    101 102 108;
    102 109 108;
    102 103 109;
    103 110 109;
    103 104 110;
    104 111 110;
    104 105 111;
    105 112 111;
    105 106 112;
    106 113 112;
    106 15 113;
    15 16 113;
    %8
    58 107 57;
    107 114 57;
    107 108 114;
    108 115 114;
    108 109 115;
    109 116 115;
    109 110 116;
    110 117 116;
    110 111 117;
    111 118 117;
    111 112 118;
    112 119 118;
    112 113 119;
    113 120 119;
    113 16 120;
    16 17 120;
    %9
    49 50 48;
    50 121 48;
    50 51 121;
    51 122 121;
    51 52 122;
    52 123 122;
    52 53 123;
    53 124 123;
    53 54 124;
    54 125 124;
    54 55 125;
    55 126 125;
    55 56 126;
    56 127 126;
    56 57 127;
    57 128 127;
    57 114 128;
    114 129 128;
    114 115 129;
    115 130 129;
    115 116 130;
    116 131 130;
    116 117 131;
    117 132 131;
    117 118 132;
    118 133 132;
    118 119 133;
    119 134 133;
    119 120 134;
    120 135 134;
    120 17 135;
    17 18 135;
    %10
    48 121 47;
    121 136 47;
    121 122 136;
    122 137 136;
    122 123 137;
    123 138 137;
    123 124 138;
    124 139 138;
    124 125 139;
    125 140 139;
    125 126 140;
    126 141 140;
    126 127 141;
    127 142 141;
    127 128 142;
    128 143 142;
    128 129 143;
    129 144 143;
    129 130 144;
    130 145 144;
    130 131 145;
    131 146 145;
    131 132 146;
    132 147 146;
    132 133 147;
    133 148 147;
    133 134 148;
    134 149 148;
    134 135 149;
    135 150 149;
    135 18 150;
    18 19 150;
    %11
    47 136 46;
    136 151 46;
    136 137 151;
    137 152 151;
    137 138 152;
    138 153 152;
    138 139 153;
    139 154 153;
    139 140 154;
    140 155 154;
    140 141 155;
    141 156 155;
    141 142 156;
    142 157 156;
    142 143 157;
    143 158 157;
    143 144 158;
    144 159 158;
    144 145 159;
    145 160 159;
    145 146 160;
    146 161 160;
    146 147 161;
    147 162 161;
    147 148 162;
    148 163 162;
    148 149 163;
    149 164 163;
    149 150 164;
    150 165 164;
    150 19 165;
    19 20 165;
    %12
    46 151 45;
    151 166 45;
    151 152 166;
    152 167 166;
    152 153 167;
    153 168 167;
    153 154 168;
    154 169 168;
    154 155 169;
    155 170 169;
    155 156 170;
    156 171 170;
    156 157 171;
    157 172 171;
    157 158 172;
    158 173 172;
    158 159 173;
    159 174 173;
    159 160 174;
    160 175 174;
    160 161 175;
    161 176 175;
    161 162 176;
    162 177 176;
    162 163 177;
    163 178 177;
    163 164 178;
    164 179 178;
    164 165 179;
    165 180 179;
    165 20 180;
    20 21 180;
    %13
    45 166 44;
    166 181 44;
    166 167 181;
    167 182 181;
    167 168 182;
    168 183 182;
    168 169 183;
    169 184 183;
    169 170 184;
    170 185 184;
    170 171 185;
    171 186 185;
    171 172 186;
    172 187 186;
    172 173 187;
    173 188 187;
    173 174 188;
    174 189 188;
    174 175 189;
    175 190 189;
    175 176 190;
    176 191 190;
    176 177 191;
    177 192 191;
    177 178 192;
    178 193 192;
    178 179 193;
    179 194 193;
    179 180 194;
    180 195 194;
    180 21 195;
    21 22 195;
    %14
    44 181 43;
    181 196 43;
    181 182 196;
    182 197 196;
    182 183 197;
    183 198 197;
    183 184 198;
    184 199 198;
    184 185 199;
    185 200 199;
    185 186 200;
    186 201 200;
    186 187 201;
    187 202 201;
    187 188 202;
    188 203 202;
    188 189 203;
    189 204 203;
    189 190 204;
    190 205 204;
    190 191 205;
    191 206 205;
    191 192 206;
    192 207 206;
    192 193 207;
    193 208 207;
    193 194 208;
    194 209 208;
    194 195 209;
    195 210 209;
    195 22 210;
    22 23 210;
    %15
    43 196 42;
    196 211 42;
    196 197 211;
    197 212 211;
    197 198 212;
    198 213 212;
    198 199 213;
    199 214 213;
    199 200 214;
    200 215 214;
    200 201 215;
    201 216 215;
    201 202 216;
    202 217 216;
    202 203 217;
    203 218 217;
    203 204 218;
    204 219 218;
    204 205 219;
    205 220 219;
    205 206 220;
    206 221 220;
    206 207 221;
    207 222 221;
    207 208 222;
    208 223 222;
    208 209 223;
    209 224 223;
    209 210 224;
    210 225 224;
    210 23 225;
    23 24 225;
    %16
    42 211 41;
    211 40 41;
    211 212 40;
    212 39 40;
    212 213 39;
    213 38 39;
    213 214 38;
    214 37 38;
    214 215 37;
    215 36 37;
    215 216 36;
    216 35 36;
    216 217 35;
    217 34 35;
    217 218 34;
    218 33 34;
    218 219 33;
    219 32 33;
    219 220 32;
    220 31 32;
    220 221 31;
    221 30 31;
    221 222 30;
    222 29 30;
    222 223 29;
    223 28 29;
    223 224 28;
    224 27 28;
    224 225 27;
    225 26 27;
    225 24 26;
    24 25 26];
u_output=u;
% u_output=[u_output; u(1,1); u(11,1); u(21,1)];
% u_output=[u_output; u(N_x+1,1); u(N_x+11,1); u(N_x+21,1)];
% u_output=[u_output; u(2*N_x+1,1); u(2*N_x+11,1); u(2*N_x+21,1)];
% u_output=[u_output; u(3*N_x+1,1); u(3*N_x+11,1); u(3*N_x+21,1)];

for i=N_nodes+1:N_nodes+161
    xp=vortices(i,1);
    yp=vortices(i,2);
    for j=1:N_nodes
        beta=x_elements(j+1,1:2)-x_elements(j,1:2);
        alpha=x_elements(j,1:2)-[xp yp];
        %g vector
        for k=1:N_p
            tmp2(k)=log(norm(alpha+beta*p_new(k)));
        end
        g_vector(j)=-(L(j)/2/pi)*(w_new*tmp2');
            
        %H matrix
        for k=1:N_p
            tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
        end
        h_vector(j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
    end
    u_output=[u_output; g_vector*q-h_vector*u];
end

patch('Vertices',vortices,'Faces',B,'FaceVertexCData',u_output,'FaceColor','interp');
    
        
        
        
    





    


        
        
