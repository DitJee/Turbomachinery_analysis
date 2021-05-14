% Written by Dit Dejphachon
% This code calculates the important paramters used in designing compressor
function [alpha_2,alpha_3,Beta_1,Beta_2,T_t3,P_t3,M_3,T_3,P_3,As] = compressorCal(T_t1,P_t1,omega,r,alpha_1,alpha_3,mdot,M_1,M_3,u_2_over_u_1,P_t3_over_P_t1)
% %% fluid's constants
% gamma = 1.4;
% c_p = 0.24; %Btu/(lbm*degreeR)
% R_g_c = 1716; %ft^2/(s^2*degree R)
% c_pg_c = 6006; %ft^2/(s^2*degree R)
% R = 151; % ft*lbf/(lbm*degree R)
% g_c =32.174;
%% fluid's constants
gamma = 1.4;
c_p = 1004; %Btu/(lbm*degreeR)
R_g_c = 287; %ft^2/(s^2*degree R)
c_pg_c = 1004; %ft^2/(s^2*degree R)
R = 287; % ft*lbf/(lbm*degree R)
g_c =1;
%% assuming rotor -> stator -> stator

%% step 1 (rotor)

T_1 = T_t1/(1+((gamma-1)/2)*M_1^2);
a_1 = sqrt(gamma*R_g_c*T_1);
V_1 = M_1*a_1;
u_1 = V_1*cos(deg2rad(alpha_1));
v_1 = V_1*sin(deg2rad(alpha_1));
P_1 = P_t1/(1+((gamma-1)/2)*M_1^2)^(gamma/(gamma-1));

MFP_1 = sqrt(gamma*g_c/R)*M_1*(1+((gamma-1)/2)*M_1^2)^((gamma+1)/(2*(gamma-1)));
A_1 = (mdot*sqrt(T_t1))/(P_t1*cos(deg2rad(alpha_1))*MFP_1);
%omegar = omega*(r/r);
omegar = omega*(r);
v_1R = omegar-v_1;
Beta_1 = rad2deg(tan(v_1R/u_1)^(-1));
V_1R = sqrt(u_1^2+v_1R^2);
M_1R = V_1R/a_1;
T_t1R = T_1*(1+((gamma-1)/2)*M_1R^2);
P_t1R = P_1*(T_t1R/T_1)^(gamma/(gamma-1));
P_t2 = P_t3_over_P_t1*P_t1;
T_t2 = T_t1*(P_t2/P_t1)^((gamma-1)/gamma);
T_t3 = T_t2;
tanBeta_2 = ((u_2_over_u_1)^(-1))*(tan(deg2rad(Beta_1))-((c_pg_c/(omegar*u_1))*(T_t2-T_t1)));
Beta_2 = rad2deg(atan(tanBeta_2));
u_2 = u_2_over_u_1*u_1;
v_2R = u_2*tanBeta_2;
V_2R = sqrt((u_2^2)+(v_2R^2));
v_2 = omegar - v_2R;
alpha_2 = rad2deg(atan(v_2/u_2));
V_2 = sqrt((u_2^2)+(v_2^2));
T_t2R = T_t1R;
P_t2R = P_t1R;

%% step 2 (stator#1)

T_2 = T_t2-((V_2^2)/(2*c_pg_c));
P_2 = P_t2*(T_2/T_t2)^(gamma/(gamma-1));
a_2 = sqrt(gamma*R_g_c*T_2);
M_2 = V_2/a_2;
M_2R = V_2R/a_2;
MFP_2 = sqrt(gamma*g_c/R)*M_2*(1+((gamma-1)/2)*M_2^2)^((gamma+1)/(2*(gamma-1)));
A_2 = (mdot*sqrt(T_t2))/(P_t2*cos(deg2rad(alpha_2))*MFP_2);

%% step 3 (stator#2)

T_3 = T_t3/(1+((gamma-1)/2)*M_3^2);
P_3 = (P_t3_over_P_t1*P_t1)*(T_3/T_t3)^(gamma/(gamma-1));
P_t3 = P_3/((T_3/T_t3)^(gamma/(gamma-1)));
a_3 = sqrt(gamma*R_g_c*T_3);
V_3 = M_3*a_3;
u_3 = V_3*cos(deg2rad(alpha_3));
v_3 = V_3*sin(deg2rad(alpha_3));
M_3 = V_3/a_3;
MFP_3 = sqrt(gamma*g_c/R)*M_3*(1+((gamma-1)/2)*M_3^2)^((gamma+1)/(2*(gamma-1)));
A_3 = (mdot*sqrt(T_t3))/((P_t3_over_P_t1*P_t1)*cos(deg2rad(alpha_3))*MFP_3);

%% package
As = [];
if A_3 > 0
    if A_2 < 0
%         A_2 = A_3*(52.3752/64.9075);
        A_2 = (A_3+A_1)./2;
    end
    As = [A_1 A_2 A_3];
end 


end

