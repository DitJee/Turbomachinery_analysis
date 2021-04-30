% Written by Dit Dejphachon
% This code calculates the important paramters used in designing compressor
function [Beta_2,Beta_3,M_3,P_3,P_t3,T_3,T_t3,A_2,A_3] = turbineCal(T_t1,P_t1,alpha_2,alpha_3,M_2,omegar,gamma,R,g_c,mdot)
% a stage consists of a stator and a rotor
% stator
T_t2 = T_t1;
T_2 = T_t2/(1+((gamma-1)/2)*M_2^2);
g_cc_p = g_c*((gamma)/(gamma-1))*R;
V_2 = sqrt(2*g_cc_p*(T_t2-T_2));
u_2 = V_2*cos(deg2rad(alpha_2));
v_2 = V_2*sin(deg2rad(alpha_2));
v_2R = v_2 - omegar;
V_2R = sqrt((u_2^2)+(v_2R^2));
Beta_2 = rad2deg(atan(v_2R/u_2));
M_2R = M_2*(V_2R/V_2);
T_t2R = T_2 + ((V_2R^2)/(2*g_cc_p));
%rotor
v_3 = 0;
u_3 = u_2;
V_3 = u_3;
v_3R = v_3 + omegar;
V_3R = sqrt((u_3^2)+(v_3R^2));
Beta_3 = rad2deg(atan(v_3R/u_3));
T_t3 = T_t2-((omegar/g_cc_p)*(v_3+v_2));
T_3 = T_t3-((V_3^2)/(2*g_cc_p));
M_3 = sqrt((2/(gamma-1))*((T_t3/T_3)-1));
M_3R = M_3*(V_3R/V_3);
T_t3R = T_t2R;
P_t2 = P_t1;
P_2 = P_t2*(T_2/T_t2)^(gamma/(gamma-1));
P_t2R = P_2*(T_t2R/T_2)^(gamma/(gamma-1));
P_t3R = P_t2R;
P_t3 = P_t2*(T_t3/T_t2)^(gamma/(gamma-1));
P_3 = P_t3*(T_3/T_t3)^(gamma/(gamma-1));
% station area
MFP_2 = sqrt(gamma*g_c/R)*M_2*(1+((gamma-1)/2)*M_2^2)^((gamma+1)/(2*(gamma-1)));
A_2 = (mdot*sqrt(T_t2))/(P_t2*cos(deg2rad(alpha_2))*MFP_2);
MFP_3 = sqrt(gamma*g_c/R)*M_3*(1+((gamma-1)/2)*M_3^2)^((gamma+1)/(2*(gamma-1)));
A_3 = (mdot*sqrt(T_t3))/(P_t3*cos(deg2rad(alpha_3))*MFP_3);
end

