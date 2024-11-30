clear all;
alpha0=0;
%  [Pos,Cl_dis]=read_lift_dis("lift_interval_10");
 Cl_dis = [2.723027e+01, 2.722310e+01, 2.720238e+01, 2.716807e+01, 2.711913e+01, 2.705459e+01, 2.697266e+01, 2.687125e+01, 2.674755e+01, 2.659755e+01, 2.641710e+01, 2.619901e+01, 2.593162e+01, 2.560213e+01, 2.519707e+01, 2.469378e+01, 2.406951e+01, 2.326630e+01, 2.219520e+01, 2.063614e+01, 1.830834e+01, 1.506318e+01]/(0.5*1.225*400*0.3);
N=length(Cl_dis);
 
 Pos=linspace(0,1.45,N);
center=(Pos(1:end-1)+Pos(2:end))/2;
len=center(2:end)-center(1:end-1);
len=[center(1),len,1.5-center(end)];
plot(Pos,Cl_dis);
hold on;
 massrate=0.12;
Ai=40/180*pi;
As=50/180*pi;
Bi=0.003;
Bs=0.007;
Gama=0;
epsilon=2/3;
vinj=massrate/1.225/Bi;
vsuc=massrate/1.225/Bs;
[Xl, Yl] = meshgrid(4:-2:-4, [0.03,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.135,0.15,0.16,0.17,0.18,0.2,0.22,0.25,0.275,0.3]); 
[Xd, Yd] = meshgrid(4:-2:-4, [0.05,0.1,0.15,0.2]); 
cl_map=table2array(readtable('cl_map_without_jet.xlsx'));
cd_map=table2array(readtable('cd_map_without_jet.xlsx'));
cdp_map=table2array(readtable('cdp_map.xlsx'));
cdf_map=table2array(readtable('cdf_map.xlsx'));
Fx_map=table2array(readtable('jetFx_map.xlsx'));
Fy_map=table2array(readtable('jetFy_map.xlsx'));
Cd_dis=zeros(N,1);
new_Cl_dis=zeros(N,1);
new_Cd_dis=zeros(N,1);
Alpha_i=zeros(N,1);
for i=1:N
    cl=Cl_dis(i)/(cos(Gama)^2);
    alpha_i=0;
    last_alpha_i=0;
    cd_eff=0;   
    while alpha_i==0 || abs(last_alpha_i-alpha_i)>0.001
        last_alpha_i=alpha_i;
        cl_eff=cl*(cos(alpha_i))+cd_eff*tan(alpha_i);
        [alpha_eff,cd_eff]=search_in_map(cl_eff);
        alpha_i=alpha0/cos(Gama)-alpha_eff;      
    end
    Alpha_i(i)=alpha_i;
    Cd_dis(i)=cd_eff/cos(alpha_i)+cl_eff*tan(alpha_i)/cos(alpha_i);
end
DeltaCl=zeros(round(epsilon*length(Cl_dis)),1);
DeltaCd=zeros(round(epsilon*length(Cl_dis)),1);
DeltaCd_p=zeros(round(epsilon*length(Cl_dis)),1);
DeltaCd_f=zeros(round(epsilon*length(Cl_dis)),1);
Cl=0;
Cd=0;
for i=1:length(Cl_dis)
    if Pos(i)<epsilon*1.5
 
        DeltaCl(i)=interp2(Xl, Yl, cl_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate, 'linear');
        DeltaCd_p(i)=interp2(Xd, Yd, cdp_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate, 'linear');
        DeltaCd_f(i)=interp2(Xd, Yd, cdf_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate, 'linear');
        DeltaCd(i)=DeltaCd_p(i)+DeltaCd_f(i);
        %垂直1/4弦长处
        new_Cl_dis(i)=DeltaCl(i)/cos(Alpha_i(i))-DeltaCd(i)*tan(Alpha_i(i))/cos(Alpha_i(i));
        new_Cd_dis(i)=DeltaCd(i)/cos(Alpha_i(i))+DeltaCl(i)*tan(Alpha_i(i))/cos(Alpha_i(i));
%         new_Cl_dis(i)=new_Cl_dis(i)-interp2(X,Y,Fy_map,(alpha0)/pi*180,massrate,'linear')/(0.5*1.225*400*0.3);
%         new_Cd_dis(i)=new_Cd_dis(i)-interp2(X,Y,Fx_map,(alpha0)/pi*180,massrate,'linear')/(0.5*1.225*400*0.3);
%         new_Cl_dis(i)=(DeltaCd_p(i)*cos(Gama)^3+DeltaCd_f(i))/cos(Alpha_i(i))+DeltaCl(i)*tan(Alpha_i(i))/cos(Alpha_i(i));
    else
        new_Cl_dis(i)=Cl_dis(i);
        new_Cd_dis(i)=Cd_dis(i);
%         Cl=Cl+Cl_dis(i);
%         Cd=Cd+Cd_dis(i);
    end

end
Cl=len*new_Cl_dis/1.5;
Cd=len*new_Cd_dis/1.5;
Cl=Cl-interp2(Xl,Yl,Fy_map,(alpha0)/pi*180,massrate,'linear')/(0.5*1.225*400*0.45);
Cd=Cd-interp2(Xl,Yl,Fx_map,(alpha0)/pi*180,massrate,'linear')/(0.5*1.225*400*0.45);
plot(Pos,new_Cl_dis);

cfd_lift= [4.980687e+01, 4.995161e+01, 5.026819e+01, 5.063544e+01, 5.082359e+01, 5.061597e+01, 5.004650e+01, 4.908667e+01, 4.764337e+01, 4.571406e+01, 4.341806e+01, 4.257368e+01, 4.219295e+01, 4.099776e+01, 3.791629e+01, 3.059722e+01, 2.846992e+01, 2.671473e+01, 2.493405e+01, 2.286124e+01, 2.033951e+01, 1.729018e+01]/(0.5*1.225*400*0.3);
plot(Pos,abs(cfd_lift));
 

% [Pos,Cl_dis]=read_lift_dis("mass_rate_0.1_alpha_0");
% sum(DeltaCl)*round(epsilon*length(Cl_dis))/length(Cl_dis)
% sum(DeltaCd)*round(epsilon*length(Cl_dis))/length(Cl_dis)
 
function [Alpha_eff,Cd_p,Cd_f]=search_in_map(Cl_eff)
    Cl_map=[0.104,0.296,0.485,0.667234041,0.837096513];
    alpha_map=[-4,-2,0]/180*pi;
    Cdp_map=[0.007292155,0.007278948,0.008104288,0.009828387,0.012613109];
    Cdf_map=[0.010385784,0.010660165,0.010845842,0.010943995,0.010955707];
    N=length(Cl_map);
    for i=1:N
        if Cl_eff<Cl_map(i)
            Cd_p=Cdp_map(i-1)+(Cl_eff-Cl_map(i-1))/(Cl_map(i)-Cl_map(i-1))*(Cdp_map(i)-Cdp_map(i-1));
            Cd_f=Cdf_map(i-1)+(Cl_eff-Cl_map(i-1))/(Cl_map(i)-Cl_map(i-1))*(Cdf_map(i)-Cdf_map(i-1));
            Alpha_eff=alpha_map(i-1)+(Cl_eff-Cl_map(i-1))/(Cl_map(i)-Cl_map(i-1))*(alpha_map(i)-alpha_map(i-1));
            break;
        end
    end
end


function [pos,local_lift]=read_lift_dis(filename)
%      filename = 'lift.txt'; % 替换为你的文件名
    
    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件');
    end
    
    % 初始化存储数据的变量
    data = [];
    keyLabel = '';
    
    % 逐行读取文件内容
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'xy/key/label') % 检测到 "xy/key/label"
            keyLabel = extractBetween(line, '"', '"');
        elseif ~isempty(regexp(line, '^\d', 'once')) % 检测数字开头的行
            % 将行内容转换为数值数组
            nums = sscanf(line, '%f');
            data = [data; nums']; % 添加到数据矩阵中
        end
    end
    
    % 关闭文件
    fclose(fid);
    
    % 显示解析的结果
    disp(['Key Label: ', keyLabel]);
    disp('Extracted Data:');
    disp(data);
    
    % 绘图
%     figure;
    % plot(data(:,1), data(:,2), '-o');
    % hold on;
     N=length(data(2:end,2))-1;
    local_lift=(data(2:end,2)-data(1:end-1,2))/(0.5*1.225*400*0.3*1.5/N); 
    plot(data(2:end,1), local_lift);
    title('Cumulative Force');
    xlabel('Distance (m)');
    ylabel('Cumulative Force (N)');
    hold on;
    pos=data(2:end,1);
end



