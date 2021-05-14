clear
close all
clc
disp(newline)
% Written by Dit Dejphachon
% This code determines the number of stages used in the compressor to get
% the desired pressure to the combustor and the properties of the
% components
%% input vars
% T_t1 = 545.61; %K
% P_t1 = 171922.9548; %pa
% omega = 8000; %rad/s
% r = 12*2.54/100; %m
% alpha_1 = 40; %degree
% alpha_3 = alpha_1; %degree
% mdot = 90; %kg/s
% M_1 = 0.3;
% M_3 = M_1;
% u_2_over_u_1 = 1.1;
% P_t3_over_P_t1 = 1.3;
%% input vars
T_t1 = 340.4; %K
P_t1 = 0.911807*100000; %pa
omega = 8000; %rad/s
r = 12*2.54/100; %m
alpha_1 = 40; %degree

alpha_3 = alpha_1; %degree
mdot = 30; %kg/s
M_1 = 0.3;
M_3 = M_1;
u_2_over_u_1 = 1.1;
P_t3_over_P_t1 = 1.3;
%% conditions
P_wanted = P_t1*9;
N_stage = 1;
Compressor_params = [];
Areas = [];
report = [];
%% compressor analysis
% assume isentropic flow
while 1
    [for_col_1,for_col_1R,for_col_2R,for_col_2,for_col_3,alpha_2,alpha_3,Beta_1,Beta_2,T_t3,P_t3,M_3,T_3,P_3,As] = compressorCal(T_t1,P_t1,omega,r,alpha_1,alpha_3,...
    mdot,M_1,M_3,u_2_over_u_1,P_t3_over_P_t1);
    output = [alpha_2,alpha_3,Beta_1,Beta_2,T_t3,P_t3,M_3,T_3,P_3];
    report = [report; for_col_1 for_col_1R for_col_2R for_col_2 for_col_3];
    Compressor_params = [Compressor_params; output];
    Areas = [Areas;As];
    if P_3 > P_wanted
        disp([' Total number of stage is ',num2str(N_stage), ' with final pressure of ', num2str(P_3), ' pa'])
        break
    end
    N_stage = N_stage+1;
    T_t1 = T_t3; %R
    P_t1 = P_t3; %psia
    alpha_1 = 40; %degree
    alpha_3 = alpha_1; %degree
    M_1 = M_3;
    M_3 = M_1;
    
end

%% summary the output(number of stage and its paramters)
disp(newline)
disp("Output parameters in each stage")
disp("  alpha_2    alpha_3    Beta_1    Beta_2    T_t3     P_t3       M_3      T_3        P_3")
disp(Compressor_params)
disp("  A_1      A_2      A_3")
disp(Areas)

%% radius calculation
% n stages before exit
rs = [];
r_h_over_r_t = 0.65;
r_t = sqrt(Areas(1:(size(Areas,1)-1),:)./(pi*(1-(1-r_h_over_r_t)^2)));
r_h = r_t.*r_h_over_r_t;
r_m = (r_t+r_h)./2;
rs = [rs;r_t r_h r_m];
% stage before exit
r_h_over_r_t = (0.92+0.9)/2;
r_t = sqrt(Areas((size(Areas,1)),:)./(pi*(1-(1-r_h_over_r_t)^2)));
r_h = r_t.*r_h_over_r_t;
r_m = (r_t+r_h)./2;
rs = [rs;r_t r_h r_m];
disp(newline)
disp("Output radius in each stage")
disp("  A_1_r_t    A_2_r_t    A_3_r_t   A_1_r_h   A_2_r_h   A_3_r_h  A_1_r_m   A_2_r_m   A_3_r_m")
disp(rs)

%% compressor height calculation
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
%c_over_h = 0.9;
c_over_h = 0.9;
W_r = ((hs(:,2)+hs(:,1))./2).*(c_over_h).*cos(deg2rad(alpha_1));
W_s = ((hs(:,3)+hs(:,2))./2).*(c_over_h).*cos(deg2rad(Compressor_params(:,2)));
Ls = [Ls;W_r W_s]
%Lenght_total = N_stage.*((sum(W_r+(1.*W_s)))+(sum(W_r+(1.*W_s)))./4)
Lenght_total = ((sum(W_r+(1.*W_s)))+(sum(W_r+(1.*W_s)))./4)
%plot
% figure()
% plot(c_over_h,Lenght_total)
% xlabel('c/h')
% ylabel('Lenght_total')

%% blade properties
% the value of c_over_h must be optimize before head
%c_over_h = 0.15;
blade_c = c_over_h.*hs;
blade_s = 1./blade_c;
temp_1 = (2.*pi.*rs(:,7:9));
number_of_blade = ceil(temp_1.*blade_s);
blade_properties = [blade_c number_of_blade];
disp(newline)
disp("Output blade properties in each stage")
disp("   b_c_A_1    b_c_A_2   b_c_A_3  numb_A_1  numb_A_2  numb_A_3")
disp(blade_properties)

%% for report
disp(newline)
disp("report")
disp("   b_c_A_1    b_c_A_2   b_c_A_3  numb_A_1  numb_A_2  numb_A_3")
disp(report)

