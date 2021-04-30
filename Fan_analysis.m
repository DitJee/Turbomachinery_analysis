clear
close all
clc
disp(newline)
% Written by Dit Dejphachon
% This code determines the number of stages used in the fan to get
% the desired pressure to the compressor and the properties of the
% components

%% input vars
T_t1 = 1345.392; %R
P_t1 = 156; %psia
omega = 1000; %rad/s
r = 12; %in
alpha_1 = 40; %degree
alpha_3 = alpha_1; %degree
mdot = 50; %lbm/s
M_1 = 0.86;
M_3 = M_1;
u_2_over_u_1 = 1.1;
P_t3_over_P_t1 = 1/3;
%% conditions
P_wanted = P_t3_over_P_t1*P_t1;
N_stage = 1;
Compressor_params = [];

%% fan analysis
% assume isentropic flow
while 1
    [alpha_2,Beta_1,Beta_2,T_t2,P_t2,M_2,T_2,P_2,A_1,A_2] = fanCal(T_t1,P_t1,omega,r,alpha_1,alpha_3,...
        mdot,M_1,M_3,u_2_over_u_1,P_t3_over_P_t1);
    output = [alpha_2,Beta_1,Beta_2,T_t2,P_t2,M_2,T_2,P_2,A_1,A_2];
    Compressor_params = [Compressor_params; output];
    if P_2 < P_wanted
        disp([' Total number of stage is ',num2str(N_stage), ' with final pressure of ', num2str(P_2), ' psia'])
        break
    end
    N_stage = N_stage+1;
    T_t1 = T_t2; %R
    P_t1 = P_t2; %psia
    r = 12; %in
    alpha_1 = 40; %degree
    alpha_3 = alpha_1; %degree
    mdot = 50; %lbm/s
    M_1 = M_2;
    M_3 = M_1;
   
    
end

%% summary the output(number of stage and its paramters)
disp(newline)
disp("Output parameters in each stage")
disp("  alpha_2    Beta_1    Beta_2   T_t2       P_t2        M_2     T_2        P_2     A_1       A_2")
disp(Compressor_params)

%% area in turbine

Areas = [Compressor_params(:,9:10)]

%% turbine radius calculation
rs = [];
% r_h_over_r_t must be adjusted 
r_h_over_r_t = 0.2;
r_t = sqrt(Areas./(pi.*(1-(1-r_h_over_r_t).^2)));
r_h = r_t.*r_h_over_r_t;
r_m = (r_t+r_h)./2;
rs = [rs;r_t r_h r_m];
disp(newline)
disp("Output radius in each stage")
disp("  A_1_r_t    A_2_r_t   A_1_r_h   A_2_r_h  A_1_r_m   A_2_r_m")
disp(rs)

%% turbine height calculation
hs = [];
h_A_1 = rs(:,1)-rs(:,3);
h_A_2 = rs(:,1+1)-rs(:,3+1);
hs = [hs;h_A_1 h_A_2];
disp(newline)
disp("Output height in each stage")
disp("    h_A_1     h_A_2")
disp(hs)

%% compressor length calculation
Ls = [];
c_over_h = 0.5;
W_r = ((hs(:,2)+hs(:,1))./2).*(c_over_h).*cos(deg2rad(alpha_1));
Lenght_total = N_stage.*((W_r))

%% blade properties
% the value of c_over_h must be optimize before head
%c_over_h = 0.8;
blade_c = c_over_h.*hs;
blade_s = 1./blade_c;
temp_1 = (2.*pi.*rs(:,5:6));
number_of_blade = ceil(temp_1.*blade_s);
blade_properties = [blade_c number_of_blade];
disp(newline)
disp("Output blade properties in each stage")
disp("   b_c_A_1    b_c_A_2  numb_A_1  numb_A_2")
disp(blade_properties)