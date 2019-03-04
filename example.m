%% ----------------------------------------
%clearing variables and closing all windows
clear
close all
%% -----------------------------------
%known variables are given their value
% Aileron properties
%thet = [0,26,-26];
%for the = 1:length(thet)
c_a = 0.547;                    %[m]
l_a = 2.771;                    %[m]
x_1 = 0.153;                    %[m]
x_2 = 1.281;                    %[m]
x_3 = 2.681;                    %[m]
x_a = 0.280;                    %[m]
x_Plocked = x_2 - (x_a/2);      %[m]
x_Pactuator = x_2 + (x_a/2);    %[m]
h_a = 0.225;                    %[m]
r_a = h_a/2;                    %[m]
t_sk = 0.0011;                  %[m]
t_sp = 0.0029;                  %[m]
t_st = 0.0012;                  %[m]
h_st = 0.015;                   %[m]
w_st = 0.020;                   %[m]
n_st = 17;                      %[-]
d_1 = 0.01103;                  %[m]    
d_3 = 0.01642;                  %[m]
theta = 0;                      %[deg]
P_actuator = 91.7*10^3;         %[N]
q = 4.53*10^3;                  %[N/m]  

% Material properties
E = 73.1*10^9;      %[Pa]
g = 28*10^9;        %[Pa]
v = 0.33;           %[-]
tau_u = 283*10^6;   %[Pa]
tau_v = 283*10^6;   %[Pa]
rho = 2.78;         %[kg/m^3]

%% -------------------------------------------------------
%in the following part of the code the different functions
%the function are diveded in such a way that their results are meaningful
%and useful for other calculations

%% Geometric Properties
% Cross Section Geometry
cross_section_geometry = @CSG;
[a_stringer,l_skin1,z_spar,l_le,l_te,l_spacing] = cross_section_geometry(t_st,w_st,h_st,c_a,h_a,n_st);

% Location Booms
location_booms = @LB;
[y_b,z_b] = location_booms(l_spacing,l_le,r_a,z_spar,l_te,c_a);

% Centroid Location
centroid_location = @CL;
[z_centroid,y_centroid,z_skin1,z_skin2,z_spar,a_spar,a_skin1,a_skin2] = centroid_location(c_a,r_a,t_sk,l_skin1,t_sp,h_a,a_stringer,z_b,n_st);

% Moment of Inertia (analytical)
moment_of_inertia = @MOI;
[I_zz,I_yy] = moment_of_inertia(t_sk,r_a,t_sp,h_a,a_skin1,a_stringer,a_skin2,z_skin2,z_centroid,a_spar,z_skin1,y_b,z_spar,z_b);

% Boom Area
boom_area = @BA;
[b_y,b_z,y_b_c,z_b_c] = boom_area(r_a,z_spar,y_b,z_b,y_centroid,z_centroid,a_stringer,t_sk,l_spacing,t_sp,h_a);

%% Beam Theory
beam_theory = @BT;
[Fn, m_y_,m_z_,shear_y,shear_z,dislist,v_y,v_z] = beam_theory(x_1,x_2,x_3,x_a,h_a,theta,P_actuator,q,l_a,c_a,E,I_zz,I_yy,d_1,d_3, x_Plocked, x_Pactuator);

%% Shear
% Shear Center
shear_center = @SC;
[m_spar_y_tot,am_1,am_2,sc_z_c,q_b_y,m_c_y_tot] = shear_center(b_y,shear_y,I_zz,y_b_c,l_spacing,z_b_c,h_a,r_a,z_spar,z_centroid,l_le,t_sk,l_te);

% Shearflow Y
shearflow_y = @SFY;
[sf_y,q_t_y] = shearflow_y(m_spar_y_tot,am_1,am_2,shear_y,sc_z_c,q_b_y,g,l_le,l_spacing,l_te,h_a,t_sk);

% Shearflow Z
shearflow_z = @SFZ;
[sf_z,q_t_z,m_spar_z_tot,m_c_z_tot] = shearflow_z(am_1,am_2,shear_z,I_yy,b_z,z_b_c,g,l_le,l_spacing,l_te,h_a,t_sk);

% Shearflow Totoal
shearflow_total = @SFT;
[q_t,q_t_max,q_t_min,dtheta_dz,tau_max,tau_min,T] = shearflow_total(q_t_y,q_t_z,sf_y,sf_z,t_sk,m_spar_z_tot,m_spar_y_tot);

theta_x = cumtrapz(dislist,dtheta_dz);
theta_tx = theta+theta_x;
m_spar_tot = m_spar_y_tot + m_spar_z_tot;
m_c_tot = m_c_y_tot + m_c_z_tot;

%Total Deflection
total_deflection = @TD;
[v_y_le,v_y_te,v_z_le,v_z_te] = total_deflection(v_y,v_z,theta_tx,c_a);
%% ----------
% Useful information is printed on the command window in the following section

% Forces
fprintf( 'tau_min ='), disp(tau_min);

shear_flow_ribs = @SFR;
[sf_r] = shear_flow_ribs(shear_z,shear_y,c_a,h_a,sc_z_c,z_b_c);

sf_ribs=[];
alpha = atan((r_a)/(z_spar));
for i = 1:length(shear_z)
syms q_13 q_21 q_32
e_1 = 0 == shear_z(i) + (z_spar)*(q_13*cos(alpha) - q_21*cos(alpha));
e_2 = 0 == shear_y(i) + (r_a)*(q_13*sin(alpha) + q_21*sin(alpha)) - q_32*h_a;
e_3 = 0 == -shear_y(i)*(z_centroid-sc_z_c(1)) + q_32*0.5*pi()*(r_a)^2;
[M,N] = equationsToMatrix([e_1, e_2, e_3],[q_13,q_21,q_32]);
sf_ribs(:,i) = vpa(linsolve(M,N));
end


%theta_z = asin(v_z/);

plot(dislist,v_y-theta_tx);
grid on
ylabel('deflection [m]')
hold on
plot(dislist,v_y);
xlabel('x [m]')
legend({'Deflection + Twist','Deflection'},'Location','east')
%end
% end
% plot(dislist,T_all,'-o');
% title('Torsion')
% xlabel('x [m]')
% ylabel('T [N]')
% legend({'Neutral','Rotation up','Rotation down'},'Location','east')
%plot(dislist,m_spar_y_tot);
% plot(dislist,m_y_);
% title('Bending moment about y-axis')
% xlabel('x [m]')
% ylabel('M_y [N]')
% legend({'Neutral'},'Location','east')
% 
% plot(dislist,v_y_);
% title('Shear in x-y plane')
% xlabel('x [m]')
% ylabel('V_y [N]')
% legend({'Neutral'},'Location','east')
% 
% plot(dislist,v_z_);
% title('Shear in x-z plane')
% xlabel('x [m]')
% ylabel('V_z [N]')
% legend({'Neutral'},'Location','east')
%plot(dislist,dtheta_dz);
%plot(dislist,theta_tx);

%% -------------------
%the following part of the code defines all the different functions
%how these functions work will be discribed within the function

% Functions
%% Geometric Properties

% Cross Section Geometry
function [a_stringer,l_skin1,z_spar,l_le,l_te,l_spacing] = CSG(t_st,w_st,h_st,c_a,h_a,n_st)
r_a  = h_a/2; %radius leading edge curvature
a_stringer = t_st*w_st+h_st*t_st; %area stringer
l_skin1 = sqrt((c_a)^2-2*c_a*r_a+2*(r_a)^2);


z_spar = c_a-r_a; %distance between trailing edge and spar 
l_le = pi()*r_a; %circumference of leading edge from the top to the bottom of the spar
l_te = sqrt(z_spar^2+r_a^2); %circumference of trailing edge from the top to the bottom of the spar

l_total = l_le+2*l_te; %circumference length of aileron
l_spacing = l_total/(n_st); %spacing between booms
end

% Location Booms
function [y_b,z_b] = LB(l_spacing,l_le,r_a,z_spar,l_te,c_a)
angle_b_2 = l_spacing/l_le*pi(); %angle between first boom and second boom
angle_b_9 = r_a/z_spar*90; %angle between z-axis and trailing edge
y_b_2 = sin(angle_b_2)*r_a; %y-location of second boom
z_b_2 = cos(angle_b_2)*r_a+z_spar; %z-location of second boom
y_b_3 = sin(angle_b_2*2)*r_a; 
z_b_3 = cos(angle_b_2*2)*r_a+z_spar;
y_b_9 = r_a*l_spacing/2/l_te;
z_b_9 = z_spar*l_spacing/2/l_te;
y_b_8 = 3*y_b_9;
z_b_8 = 3*z_b_9;
y_b_7 = 5*y_b_9;
z_b_7 = 5*z_b_9;
y_b_6 = 7*y_b_9;
z_b_6 = 7*z_b_9;
y_b_5 = 9*y_b_9;
z_b_5 = 9*z_b_9;
y_b_4 = 11*y_b_9;
z_b_4 = 11*z_b_9;
y_b_a = [y_b_2,y_b_3,y_b_4,y_b_5,y_b_6,y_b_7,y_b_8,y_b_9]; %y-location of booms with positive y
z_b_a = [z_b_2,z_b_3,z_b_4,z_b_5,z_b_6,z_b_7,z_b_8,z_b_9]; %z-location of booms with positive y
y_b_b = -y_b_a; %y-location of booms with negative y
z_b_b = z_b_a; %z-location of booms with negative y
y_b = [0,y_b_a,fliplr(y_b_b)]; %y-location of all booms
z_b = [c_a,z_b_a,fliplr(z_b_b)]; %z-location of all booms
%scatter(z_b,y_b)
%hold on;
end

% Centroid Location
function [z_centroid,y_centroid,z_skin1,z_skin2,z_spar,a_spar,a_skin1,a_skin2] = CL(c_a,r_a,t_sk,l_skin1,t_sp,h_a,a_stringer,z_b,n_st)
z_skin1 = 0.5*(c_a-r_a);
y_skin1 = 0.5*r_a;
a_skin1 = t_sk*l_skin1;
z_skin2 = ((4*r_a)/(3*pi))+c_a-r_a;
y_skin2 = 0;
a_skin2 = pi*r_a*t_sk;
z_spar = c_a-r_a;
y_spar = 0;
a_spar = t_sp*h_a;

z_centroid = (2*a_skin1*z_skin1+a_skin2*z_skin2+a_spar*z_spar...
    +sum(z_b*a_stringer))/(a_stringer*n_st+a_spar+a_skin2+2*a_skin1);
y_centroid = 0;
%scatter(z_centroid,y_centroid);
%hold on;
end

% Moment of Inertia (analytical)
% Moment of Inertia (analytical)
function [I_zz,I_yy] = MOI(t_sk,r_a,t_sp,h_a,a_skin1,a_stringer,a_skin2,z_skin2,z_centroid,a_spar,z_skin1,y_b,z_spar,z_b);
I_zzskin2 = 0.5*t_sk*(r_a)^3; %Inertia of semi circle about z
I_zzspar = t_sp*(h_a)^3/12; %Inertia of spar about z
I_zzskin1 = 2*(1/12*t_sk*(r_a)^3+a_skin1*(0.5*r_a)^2); %Inertia of slopes about z
I_zzstringer = sum((y_b).^2*a_stringer); %Inertia of stringers about z
I_zz = I_zzskin2+I_zzspar+I_zzskin1+I_zzstringer; %Sum of all I_zz
%Izz = 1.252*10^(-5)    %for testing only!

I_yyskin2 = I_zzskin2 - a_skin2*((4*r_a)/(3*pi))^2+a_skin2...
    *(z_skin2-z_centroid)^2; %Inertia of semi circle about y
I_yyspar = a_spar*(z_spar-z_centroid)^2; %Inertia of spar about y
I_yyskin1 = 2*(1/12*t_sk*(z_spar)^3+a_skin1*(z_skin1-z_centroid)^2); %Inertia of slopes about y
I_yystringer = sum((z_b-z_centroid).^2*a_stringer);
I_yy = I_yyskin2+I_yyspar+I_yyskin1+I_yystringer;

%I_yy = 9.934*10^(-5)   %for testing only!
end

% Boom Area
function [b_y,b_z,y_b_c,z_b_c] = BA(r_a,z_spar,y_b,z_b,y_centroid,z_centroid,a_stringer,t_sk,l_spacing,t_sp,h_a);
y_b_spar_1 = r_a; %y-location of boom for spar for positive y
z_b_spar_1 = z_spar; %z-location of boom for spar for positive y
y_b_spar_2 = -r_a; %y-location of boom for spar for negative y
z_b_spar_2 = z_spar; %z-location of boom for spar for negative y
%scatter([z_b_spar_1 z_b_spar_2],[y_b_spar_1 y_b_spar_2]);
y_b_c_spar_1 = y_b_spar_1-y_centroid; %y-position relative to y-centroid
z_b_c_spar_1 = z_b_spar_1-z_centroid; %z-position relative to z-centroid
y_b_c_spar_2 = y_b_spar_2-y_centroid;
z_b_c_spar_2 = z_b_spar_2-z_centroid;
y_b_c = y_b-y_centroid;
z_b_c = z_b-z_centroid;

for i = 1:length(y_b) %boom area calculations as a result of moment in y for all booms
    if i == 1 %y-position of boom 1 is zero 
    b_y(i) = a_stringer +4*t_sk*l_spacing/6; %%% NAAR KIJKEN
    elseif i == length(y_b) %boom 17 is compared to boom 16 and boom 1
    b_y(i) = (t_sk*l_spacing/6)*(4+y_b_c(1)/y_b_c(i)+y_b_c(i-1)/y_b_c(i))+ a_stringer;
    else %booms are compared to booms next to it
    b_y(i) = (t_sk*l_spacing/6)*(4+y_b_c(i+1)/y_b_c(i)+y_b_c(i-1)/y_b_c(i))+ a_stringer;
    end
end
for i = 1:length(z_b) %boom area calculations as a result of moment in z for all booms
    if i == 1 %boom 1 is compared to boom 17 and boom 2
    b_z(i) = (t_sk*l_spacing/6)*(4+z_b_c(i+1)/z_b_c(i)+z_b_c(17)/z_b_c(i))+ a_stringer;
    elseif i == length(z_b) %boom 17 is compared to boom 16 and boom 1
    b_z(i) = (t_sk*l_spacing/6)*(4+z_b_c(1)/z_b_c(i)+z_b_c(i-1)/z_b_c(i))+ a_stringer;
    else %booms are compared to booms next to it
    b_z(i) = (t_sk*l_spacing/6)*(4+z_b_c(i+1)/z_b_c(i)+z_b_c(i-1)/z_b_c(i))+ a_stringer;
    end
end

b_y(length(y_b)+1) = (t_sp*h_a/6)*(2+ y_b_spar_2/y_b_spar_1); %boom spar top is compared to boom spar bottom
b_z(length(z_b)+1) = (t_sp*h_a/6)*(2+ z_b_spar_2/z_b_spar_1); 
b_y(length(y_b)+2) = (t_sp*h_a/6)*(2+ y_b_spar_1/y_b_spar_2); %boom spar bottom is compared to boom spar top
b_z(length(z_b)+2) = (t_sp*h_a/6)*(2+ z_b_spar_1/z_b_spar_2);

y_b_c = [y_b_c, y_b_c_spar_1, y_b_c_spar_2]; %all y-positions booms in same array
z_b_c = [z_b_c, z_b_c_spar_1, z_b_c_spar_2]; %all z-positions booms in same array
end
%% Beam Theory

function [Fn, m_y_,m_z_,shear_y,shear_z,dislist,v_y,v_z] = BT(x_1,x_2,x_3,x_a,h_a,theta,P_actuator,q,l_a,c_a,E,I_zz,I_yy,d_1,d_3, x_Plocked, x_Pactuator)
Amatrix = [[1,1,1,0,0,0,0,0,0,0,0], 
          [0,0,0,1,1,1,1,0,0,0,0], 
          [0,0,0,0,-(x_2-x_1),-(x_3-x_1),-(x_2-x_1-(x_a/2)),0,0,0,0], 
          [0,0,0,0,0,0,((h_a/2)-(h_a/2)*tan(theta))*cos(theta),0,0,0,0],
          [0,(x_2-x_1),(x_3-x_1),0,0,0,0,0,0,0,0],
          [0, 0, 0, 0, 0, 0, 0, x_1, 1,0,0],
          [-(1/6)*(x_2-x_1)^3*cos(theta), 0, 0, -(1/6)*(x_2-x_1)^3*sin(theta), 0, 0, -(1/6)*((x_a/2))^3*sin(theta),x_2, 1,0,0],
          [-(1/6)*(x_3-x_1)^3*cos(theta), -(1/6)*(x_3-x_2)^3*cos(theta), 0, -(1/6)*(x_3-x_1)^3*sin(theta), -(1/6)*(x_3-x_2)^3*sin(theta), 0, -(1/6)*(x_3-x_2+(x_a/2))^3*sin(theta),x_3,1,0,0],
          [0, 0, 0, 0, 0, 0, 0, 0, 0,x_1,1],
          [-(1/6)*(x_2-x_1)^3*sin(theta), 0, 0, (1/6)*(x_2-x_1)^3*cos(theta), 0, 0, (1/6)*((x_a/2))^3*cos(theta),0, 0,x_2,1],
          [-(1/6)*(x_3-x_1)^3*sin(theta), -(1/6)*(x_3-x_2)^3*sin(theta), 0, (1/6)*(x_3-x_1)^3*cos(theta), (1/6)*(x_3-x_2)^3*cos(theta), 0, (1/6)*(x_3-x_2+(x_a/2))^3*cos(theta),0,0,x_3,1],
          ];
Cmatrix =  [[q*l_a],
           [P_actuator],
           [-P_actuator*(x_2 - x_1 + (x_a/2))],
           [P_actuator*((h_a/2) - (h_a/2)*tan(theta))*cos(theta) + q*l_a*(0.25*c_a - (h_a/2)*cos(theta))],
           [(1/2)*(l_a-x_1)^2*q-0.5*x_1^2*q],
           [-E*I_zz*(d_1*cos(theta)) - (((x_1)^4/24)*q*cos(theta))],
           [-x_2^4/24*q*cos(theta)],
           [-E*I_zz*(d_3*cos(theta)) - ((x_3)^4/24)*q*cos(theta)- P_actuator* 1/6*(x_3-x_2-(x_a/2))^3*sin(theta)],
           [-E*I_yy*(d_1*sin(theta)) - (((x_1)^4/24)*q*sin(theta))],
           [-x_2^4*q*sin(theta)],
           [-E*I_yy*(d_3*sin(theta)) - (((x_3)^4/24)*q*sin(theta)) + P_actuator* (1/6)*(x_3-x_2-(x_a/2))^3*cos(theta)]];
               

Fn = linsolve(Amatrix,Cmatrix); %Matrix A and Matrix C are solved

F_1y = Fn(1);
F_2y = Fn(2);
F_3y = Fn(3);
F_1z = Fn(4);
F_2z = Fn(5);
F_3z = Fn(6);
P_locked = Fn(7);
C_1 = Fn(8);
C_2 = Fn(9);
D_1 = Fn(10);
D_2 = Fn(11);

dx = 0.00001; %step calculation shear
%dislist = [0, dx, 2*dx, 10*dx, x_1-10*dx, x_1-dx, x_1, x_1+dx, x_Plocked-dx, x_Plocked  ,x_Plocked+dx, x_2-dx, x_2, x_2+dx, x_Pactuator-dx, x_Pactuator, x_Pactuator+dx, x_3-dx, x_3, x_3+dx, l_a-dx, l_a, l_a+dx];
dislist = [0:l_a/50:l_a]; %locations of cuts over the aileron

for i = 1:length(dislist) %moment and shear calculations for all cuts
    m_y_(i) = -(1/2)*q*(dislist(i))^2*sin(theta) + heaviside(dislist(i)-x_1)*F_1y*(dislist(i)-x_1)*sin(theta) - heaviside(dislist(i)-x_1)*F_1z*(dislist(i)-x_1)*cos(theta) - heaviside(dislist(i)-x_Plocked)*P_locked*(dislist(i)-x_Plocked)*cos(theta) + heaviside(dislist(i)-x_2)*F_2y*(dislist(i)-x_2)*sin(theta) - heaviside(dislist(i)-x_2)*F_2z*(dislist(i)-x_2)*cos(theta) + heaviside(dislist(i)-x_Pactuator)*P_actuator*(dislist(i)-x_Pactuator)*cos(theta) + heaviside(dislist(i)-x_3)*F_3y*(dislist(i)-x_3)*sin(theta) - heaviside(dislist(i)-x_3)*F_3z*(dislist(i)-x_3)*cos(theta); %moment around z for cuts
    m_z_(i) = -(1/2)*q*(dislist(i))^2*cos(theta) + heaviside(dislist(i)-x_1)*F_1y*(dislist(i)-x_1)*cos(theta) + heaviside(dislist(i)-x_1)*F_1z*(dislist(i)-x_1)*sin(theta) + heaviside(dislist(i)-x_Plocked)*P_locked*(dislist(i)-x_Plocked)*sin(theta) + heaviside(dislist(i)-x_2)*F_2y*(dislist(i)-x_2)*cos(theta) + heaviside(dislist(i)-x_2)*F_2z*(dislist(i)-x_2)*sin(theta) - heaviside(dislist(i)-x_Pactuator)*P_actuator*(dislist(i)-x_Pactuator)*sin(theta) + heaviside(dislist(i)-x_3)*F_3y*(dislist(i)-x_3)*cos(theta) + heaviside(dislist(i)-x_3)*F_3z*(dislist(i)-x_3)*sin(theta); %moment around y for cuts

    m_y_1 = -(1/2)*q*(dislist(i)+dx)^2*sin(theta) + heaviside(dislist(i)+dx-x_1)*F_1y*(dislist(i)+dx-x_1)*sin(theta) - heaviside(dislist(i)+dx-x_1)*F_1z*(dislist(i)+dx-x_1)*cos(theta) - heaviside(dislist(i)+dx-x_Plocked)*P_locked*(dislist(i)+dx-x_Plocked)*cos(theta) + heaviside(dislist(i)+dx-x_2)*F_2y*(dislist(i)+dx-x_2)*sin(theta) - heaviside(dislist(i)+dx-x_2)*F_2z*(dislist(i)+dx-x_2)*cos(theta) + heaviside(dislist(i)+dx-x_Pactuator)*P_actuator*(dislist(i)+dx-x_Pactuator)*cos(theta) + heaviside(dislist(i)+dx-x_3)*F_3y*(dislist(i)+dx-x_3)*sin(theta) - heaviside(dislist(i)+dx-x_3)*F_3z*(dislist(i)+dx-x_3)*cos(theta); %moment around z for cuts
    m_z_1 = -(1/2)*q*(dislist(i)+dx)^2*cos(theta) + heaviside(dislist(i)+dx-x_1)*F_1y*(dislist(i)+dx-x_1)*cos(theta) + heaviside(dislist(i)+dx-x_1)*F_1z*(dislist(i)+dx-x_1)*sin(theta) + heaviside(dislist(i)+dx-x_Plocked)*P_locked*(dislist(i)+dx-x_Plocked)*sin(theta) + heaviside(dislist(i)+dx-x_2)*F_2y*(dislist(i)+dx-x_2)*cos(theta) + heaviside(dislist(i)+dx-x_2)*F_2z*(dislist(i)+dx-x_2)*sin(theta) - heaviside(dislist(i)+dx-x_Pactuator)*P_actuator*(dislist(i)+dx-x_Pactuator)*sin(theta) + heaviside(dislist(i)+dx-x_3)*F_3y*(dislist(i)+dx-x_3)*cos(theta) + heaviside(dislist(i)+dx-x_3)*F_3z*(dislist(i)+dx-x_3)*sin(theta); %moment around y for cuts
    
    dm_y_(i) = m_y_(i)-m_y_1; %difference between moment around z for the cuts and close to the cuts
    dm_z_(i) = m_z_(i)-m_z_1; %difference between moment around y for the cuts and close to the cuts
    shear_z(i) = dm_y_(i)/dx; %shear for z as a result of diffirentiated moment around y
    shear_y(i) = dm_z_(i)/dx; %shear for y as a result of diffirentiated moment around z

    v_y_(i) = (1/24*dislist(i)^4*q*cos(theta) - heaviside(dislist(i)-x_1)*1/6*(dislist(i)-x_1)^3*cos(theta)*F_1y - heaviside(dislist(i)-x_1)*F_1z*1/6*(dislist(i)-x_1)^3*sin(theta) - heaviside(dislist(i)-x_1+x_a/2)*1/6*P_locked*(dislist(i)-x_1+x_a/2)^3*sin(theta) - heaviside(dislist(i)-x_2)*1/6*F_2y*(dislist(i)-x_2)^3*cos(theta) - heaviside(dislist(i)-x_2)*1/6*F_2z*(dislist(i)-x_2)^3*sin(theta) + heaviside(dislist(i)-x_2-x_a/2)*1/6*P_actuator*(dislist(i)-x_2-x_a/2)^3*sin(theta) - heaviside(dislist(i)-x_3)*1/6*F_3y*(dislist(i)-x_3)^3*cos(theta) - heaviside(dislist(i)-x_3)*1/6*F_3z*(dislist(i)-x_3)^3*sin(theta) + dislist(i)*C_1+C_2)/(-E*I_zz);
    v_z_(i) = (1/24*dislist(i)^4*q*sin(theta) - heaviside(dislist(i)-x_1)*1/6*(dislist(i)-x_1)^3*sin(theta)*F_1y + heaviside(dislist(i)-x_1)*F_1z*1/6*(dislist(i)-x_1)^3*cos(theta) + heaviside(dislist(i)-x_1+x_a/2)*1/6*P_locked*(dislist(i)-x_1+x_a/2)^3*cos(theta) - heaviside(dislist(i)-x_2)*1/6*F_2y*(dislist(i)-x_2)^3*sin(theta) + heaviside(dislist(i)-x_2)*1/6*F_2z*(dislist(i)-x_2)^3*cos(theta) - heaviside(dislist(i)-x_2-x_a/2)*1/6*P_actuator*(dislist(i)-x_2-x_a/2)^3*cos(theta) - heaviside(dislist(i)-x_3)*1/6*F_3y*(dislist(i)-x_3)^3*sin(theta) + heaviside(dislist(i)-x_3)*1/6*F_3z*(dislist(i)-x_3)^3*cos(theta) + dislist(i)*D_1+D_2)/(-E*I_yy);
    v_y(i) = -v_z_(i)*sin(theta)+v_y_(i)*cos(theta);
    v_z(i) = v_z_(i)*cos(theta)+v_y_(i)*sin(theta);
end
end
%% Shear

% Shear Center
function [m_spar_y_tot,am_1,am_2,sc_z_c,q_b_y,m_c_y_tot] = SC(b_y,shear_y,I_zz,y_b_c,l_spacing,z_b_c,h_a,r_a,z_spar,z_centroid,l_le,t_sk,l_te);
for i = 1:(length(b_y)-1) %shear flow through y for all sections between booms
    if i == 15 || i == 3 %chosen location for cuts
        q_b_y(i,:) = 0;
    else
        q_b_y(i,:) = shear_y/I_zz*b_y(i)*y_b_c(i);
    end
    m_spar_y(i,:) = l_spacing*abs(z_b_c(i)-z_b_c(19))*q_b_y(i,:); %moment calculation around boom at bottom of spar
    if i == 18
        m_c_y(i,:) = q_b_y(i,:)*h_a*z_b_c(i); %moment calculation around centroid
    else
        m_c_y(i,:) = q_b_y(i,:)*l_spacing*z_b_c(i);
    end
end
m_spar_y_tot = sum(m_spar_y); %summation moments around boom at bottom of spar
m_c_y_tot = sum(m_c_y); %summation moments around centroid

am_1 = pi()*r_a^2; %area between leading edge and spar
am_2 = r_a*z_spar; %area between trailing edge and spar

sc_z = [];
sc_z_c = [];
for i = 1:length(m_spar_y_tot) %shear center calculation for all cuts
syms q_sc1 q_sc2 etha_s

sc_1_z = 0 == ( (q_b_y(15,i)+q_sc1)*(0.5*l_le-2*l_spacing) + ((q_b_y(16,i)+q_sc1)*l_spacing) + ((q_b_y(17,i)+q_sc1)*l_spacing) + ((q_b_y(1,i)+q_sc1)*l_spacing) + (q_b_y(2,i)+q_sc1)*l_spacing + (q_b_y(3,i)+q_sc1)*(0.5*l_le-2*l_spacing) + (q_b_y(18,i)+q_sc1-q_sc2)*h_a)/(t_sk);
sc_2_z = 0 == ((l_te-5.5*l_spacing)*(q_b_y(3,i)+q_b_y(15,i)+2*q_sc2)+l_spacing*(q_b_y(4,i)+q_b_y(5,i)+q_b_y(6,i)+q_b_y(7,i)+q_b_y(8,i)+q_b_y(9,i)+q_b_y(10,i)+q_b_y(11,i)+q_b_y(12,i)+q_b_y(13,i)+q_b_y(14,i)+11*q_sc2)+(q_sc2-q_b_y(18,i)-q_sc1)*h_a)/(t_sk);

[A,B] = equationsToMatrix([sc_1_z, sc_2_z],[q_sc1,q_sc2]);
sc_z(:,i) = vpa(linsolve(A,B)); %redundant shear flows to find z-axis

c_etha = shear_y(i)*(etha_s) == 2*am_1*sc_z(1,i)+2*am_2*sc_z(2,i) + m_c_y_tot(i);
sc_z_c(i) = vpa(solve(c_etha,etha_s)); %shear center on z-axis
end
end

% Shearflow Y
function [sf_y,q_t_y,m_s] = SFY(m_spar_y_tot,am_1,am_2,shear_y,sc_z_c,q_b_y,g,l_le,l_spacing,l_te,h_a,t_sk);
sf_y = [];
for i = 1:length(m_spar_y_tot) %shear force for all cuts
syms q_s1_y q_s2_y dtheta_dz

eq_1_y = 0 == m_spar_y_tot(i) + 2*am_1*q_s1_y + 2*am_2*q_s2_y+shear_y(i)*sc_z_c(i);
eq_2_y = g*dtheta_dz == 1/(2*am_1)*((q_b_y(15,i)+q_s1_y)*(0.5*l_le-2*l_spacing) + ((q_b_y(16,i)+q_s1_y)*l_spacing) + ((q_b_y(17,i)+q_s1_y)*l_spacing) + ((q_b_y(1,i)+q_s1_y)*l_spacing) + (q_b_y(2,i)+q_s1_y)*l_spacing + (q_b_y(3,i)+q_s1_y)*(0.5*l_le-2*l_spacing) + (q_b_y(18,i)+q_s1_y-q_s2_y)*h_a)/(t_sk);
eq_3_y = g*dtheta_dz == 1/(2*am_2)*((l_te-5.5*l_spacing)*(q_b_y(3,i)+q_b_y(15,i)+2*q_s2_y)+l_spacing*(q_b_y(4,i)+q_b_y(5,i)+q_b_y(6,i)+q_b_y(7,i)+q_b_y(8,i)+q_b_y(9,i)+q_b_y(10,i)+q_b_y(11,i)+q_b_y(12,i)+q_b_y(13,i)+q_b_y(14,i)+11*q_s2_y)+(q_s2_y-q_b_y(18,i)-q_s1_y)*h_a)/(t_sk);

[C,D] = equationsToMatrix([eq_1_y, eq_2_y, eq_3_y],[q_s1_y,q_s2_y,dtheta_dz]);
sf_y(:,i) = vpa(linsolve(C,D)); %shear force in y
end

for i = 1:3
    q_t_y(i,:) = q_b_y(i+14,:)+sf_y(1,:);
end
for i = 4:6
   q_t_y(i,:) = q_b_y(i-3,:)+sf_y(1,:); 
end
for i = 7:19
   q_t_y(i,:) = q_b_y(i-3,:)+sf_y(2,:);
end
q_t_y(20,:) = q_b_y(18,:)+sf_y(1,:)+sf_y(2,:); %total shear force in y
end

% Shearflow Z
function [sf_z,q_t_z,m_spar_z_tot,m_c_z_tot] = SFZ(am_1,am_2,shear_z,I_yy,b_z,z_b_c,g,l_le,l_spacing,l_te,h_a,t_sk);
for i = 1:(length(b_z)-1) %shear flow through z for all sections between booms
    if i == 15 || i ==3 %chosen location for cuts
        q_b_z(i,:) = 0;
    else
        q_b_z(i,:) = shear_z/I_yy*b_z(i)*z_b_c(i);
    end
    m_spar_z(i,:) = l_spacing*abs(z_b_c(i)-z_b_c(19))*q_b_z(i,:); %%%% Pepijn!!!!!!!!
    if i == 18
        m_c_z(i,:) = q_b_z(i,:)*h_a*z_b_c(i); %moment calculation around centroid
    else
        m_c_z(i,:) = q_b_z(i,:)*l_spacing*(z_b_c(i));
    end 
end
m_spar_z_tot = sum(m_spar_z); %summation moments around boom at bottom of spar
m_c_z_tot = sum(m_c_z);
sf_z = [];
for i = 1:length(m_spar_z_tot) %shear force for all cuts
syms q_s1_z q_s2_z dtheta_dz

eq_1_z = 0 == m_spar_z_tot(i) + 2*am_1*q_s1_z + 2*am_2*q_s2_z;
eq_2_z = g*dtheta_dz == 1/(2*am_1)*((q_b_z(15,i)+q_s1_z)*(0.5*l_le-2*l_spacing) + ((q_b_z(16,i)+q_s1_z)*l_spacing) + ((q_b_z(17,i)+q_s1_z)*l_spacing) + ((q_b_z(1,i)+q_s1_z)*l_spacing) + (q_b_z(2,i)+q_s1_z)*l_spacing + (q_b_z(3,i)+q_s1_z)*(0.5*l_le-2*l_spacing) + (q_b_z(18,i)+q_s1_z-q_s2_z)*h_a)/(t_sk);
eq_3_z = g*dtheta_dz == 1/(2*am_2)*((l_te-5.5*l_spacing)*(q_b_z(3,i)+q_b_z(15,i)+2*q_s2_z)+l_spacing*(q_b_z(4,i)+q_b_z(5,i)+q_b_z(6,i)+q_b_z(7,i)+q_b_z(8,i)+q_b_z(9,i)+q_b_z(10,i)+q_b_z(11,i)+q_b_z(12,i)+q_b_z(13,i)+q_b_z(14,i)+11*q_s2_z)+(q_s2_z-q_b_z(18,i)-q_s1_z)*h_a)/(t_sk);

[C,D] = equationsToMatrix([eq_1_z, eq_2_z, eq_3_z],[q_s1_z,q_s2_z,dtheta_dz]);
sf_z(:,i) = vpa(linsolve(C,D)); %shear force in z
end

for i = 1:3
    q_t_z(i,:) = q_b_z(i+14,:)+sf_z(1,:);
end
for i = 4:6
   q_t_z(i,:) = q_b_z(i-3,:)+sf_z(1,:); 
end
for i = 7:19
   q_t_z(i,:) = q_b_z(i-3,:)+sf_z(2,:);
end
q_t_z(20,:) = q_b_z(18,:)+sf_z(1,:)+sf_z(2,:); %total shear force in z
end

% Shearflow Rib

function [sf_r] = SFR(shear_z,shear_y,c_a,h_a,sc_z_c,z_b_c);
syms q_13 q_21 q_32
sf_r=[];
for i = 1:length(shear_z)
alpha = atan((h_a/2)/(c_a-(h_a/2)));
eq1 = 0 == shear_z + (c_a-h_a/2)*(q_13*cos(alpha)-q_21*cos(alpha));
eq2 = 0 == shear_y + (h_a/2)*(q_13*sin(alpha)+q_21*sin(alpha)-q_32);
eq3 = 0 == q_32*(c_a*(h_a/2)+pi()*(h_a/2)^2)+shear_z*(h_a/2)+shear_y*(z_b_c(18)-sc_z_c(1));
[V,O] = equationsToMatrix([eq1, eq2, eq3],[q_13,q_21,q_32]);
sf_r(:,i) = vpa(linsolve(V,O)); %shear force in z
    
end
end

%Shearflow Total
function [q_t,q_t_max,q_t_min,dtheta_dz,tau_max,tau_min,T] = SFT(q_t_y,q_t_z,sf_y,sf_z,t_sk,m_spar_z_tot,m_spar_y_tot);
q_t = q_t_y + q_t_z;
q_t_max = max(q_t);
q_t_min = min(q_t);
dtheta_dz = sf_y(3,:) +sf_z(3,:);
tau_max = q_t_max/t_sk;
tau_min = q_t_min/t_sk;
T = -(m_spar_z_tot + m_spar_y_tot);
end

%Total Deflection
function [v_y_le,v_y_te,v_z_le,v_z_te] = TD(v_y,v_z,theta_tx,c_a);
v_y_le = v_y + 0.25*c_a*sin(-theta_tx);
v_y_te = v_y + 0.75*c_a*sin(theta_tx);
v_z_le = v_z - 0.25*c_a*cos(theta_tx);
v_z_te = v_z + 0.75*c_a*cos(theta_tx);
end