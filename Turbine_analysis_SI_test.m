clear
close all
clc
disp(newline)
% Written by Dit Dejphachon
% This code determines the number of stages used in the turbine to get
% the desired pressure to the afterburning chamber and the properties of the
% components
% Tt_4 = 1800 K
% Pt_4 =14.115 bar
% Tt_6 = 1255.81 k
% Pt_6 = 2.212 bar
%% input vars
T_t1 = 2000; %K
P_t1 = 7.878*100000; %psia
%alpha_2 = 60; %degree
alpha_2 = 9; %degree
alpha_3 = 0; %degreee
M_2 = 1.1;
omegar = 1000; %ft/s
mdot = 30; %kg/s
% u_2 = u_3
gamma = 1.3;
R = 287;% ft*lbf/(lbm*degree R)
g_c = 1;
P_t3_over_P_t1 = 1/1.3;

%% conditions
P_wanted = 1.778543*100000;
N_stage = 1;
Compressor_params = [];

%% find the flow properties for this isentropic turbine stage.
while 1
    [Beta_2,Beta_3,M_3,P_3,P_t3,T_3,T_t3,A_2,A_3] = turbineCal(T_t1,P_t1,alpha_2,...
        alpha_3,M_2,omegar,gamma,R,g_c,mdot);
    output = [Beta_2,Beta_3,M_3,P_3,P_t3,T_3,T_t3,A_2,A_3];
    Compressor_params = [Compressor_params; output];
    if P_3 < P_wanted
        disp([' Total number of stage is ',num2str(N_stage), ' with final pressure of ', num2str(P_3), ' psia'])
        break
    end
    N_stage = N_stage+1;
    T_t1 = T_t3; %R
    P_t1 = P_t3; %psia
    M_2 = M_3;
    
end

%% summary the output(number of stage and its paramters)
disp(newline)
disp("Output parameters in each stage")
disp("    Beta_2    Beta_3     M_3       P_3       P_t3      T_3      T_t3       A_2       A_3")
disp(Compressor_params)
%% area in turbine
A_1 = [0.0954;Compressor_params(1,9);Compressor_params(2,9)];
Areas = [A_1 Compressor_params(:,8:9)];
%51.8357
%% turbine radius calculation
rs = [];
%r_h_over_r_t = (0.75+0.6)/2;
r_h_over_r_t = 0.5;
r_t = sqrt(Areas./(pi*(1-(1-r_h_over_r_t)^2)));
r_h = r_t.*r_h_over_r_t;
r_m = (r_t+r_h)./2;
rs = [rs;r_t r_h r_m];
disp(newline)
disp("Output radius in each stage")
disp("  A_1_r_t    A_2_r_t    A_3_r_t   A_1_r_h   A_2_r_h   A_3_r_h  A_1_r_m   A_2_r_m   A_3_r_m")
disp(rs)

%% turbine height calculation
hs = [];
h_A_1 = rs(:,1)-rs(:,4);
h_A_2 = rs(:,1+1)-rs(:,4+1);
h_A_3 = rs(:,1+2)-rs(:,4+2);
hs = [hs;h_A_1 h_A_2 h_A_3];
disp(newline)
disp("Output height in each stage")
disp("    h_A_1     h_A_2     h_A_3 ")
disp(hs)
%% compressor length calculation
Ls = [];
c_over_h = 0.95;
W_r = ((hs(:,2)+hs(:,1))./2).*(c_over_h).*cos(deg2rad(Compressor_params(:,1)));
W_s = ((hs(:,3)+hs(:,2))./2).*(c_over_h).*cos(deg2rad(Compressor_params(:,2)));
Lenght_total = ((sum(W_r+(2.*W_s)))+(sum(W_r+(2.*W_s)))./4)

%% blade properties
% the value of c_over_h must be optimize before head

blade_c = c_over_h.*hs;
blade_s = 1./blade_c;
temp_1 = (2.*pi.*rs(:,7:9));
number_of_blade = ceil(temp_1.*blade_s);
blade_properties = [blade_c number_of_blade];
disp(newline)
disp("Output blade properties in each stage")
disp("   b_c_A_1    b_c_A_2   b_c_A_3  numb_A_1  numb_A_2  numb_A_3")
disp(blade_properties)