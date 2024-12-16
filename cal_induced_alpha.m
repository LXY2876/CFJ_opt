clear all;
allResults=[];
massrateValues = [0.05,0.1,0.12, 0.15,0.2,0.22,0.24]; % 多个质量流量值示例
induced_Angle=[];
Pos=[];
% [Pos,induced_Angle]=induced_angle_data_from_excel('induce_angle_xflr5_rectangular.csv');
N=(length(Pos)+1)/2;

for i = 1:length(massrateValues)
    currentResult = saveAeroDataToExcel(massrateValues(i),-Pos(N:-1:1)',induced_Angle(N:-1:1));
    allResults = [allResults; currentResult]; % 添加到总的结果矩阵中
end
resultTable = array2table(allResults, 'VariableNames', {'massrate', 'Cl_no_jet', 'Cd_no_jet', 'Cl', 'Cd','DeltaP'});
writetable(resultTable, 'all_results.xlsx', 'Sheet', 'Sheet1');
 
function result=saveAeroDataToExcel(real_massrate,Pos,induced_Angle)
    shapes= size(induced_Angle);
    if shapes(2)==0

    %矩阵机翼CFD
%      Cl_dis =[3.890195e+01, 3.887991e+01, 3.884007e+01, 3.878163e+01, 3.870317e+01, 3.860339e+01, 3.847709e+01, 3.832131e+01, 3.813236e+01, 3.790547e+01, 3.763276e+01, 3.730492e+01, 3.690991e+01, 3.643016e+01, 3.584208e+01, 3.509842e+01, 3.418522e+01, 3.301707e+01, 3.146745e+01, 2.927147e+01, 2.603197e+01, 2.156876e+01];
%     Cl_dis = [2.722508e+01, 2.720927e+01, 2.718091e+01, 2.713953e+01, 2.708415e+01, 2.701388e+01, 2.692510e+01, 2.681580e+01, 2.668348e+01, 2.652493e+01, 2.633479e+01, 2.610681e+01, 2.583289e+01, 2.550113e+01, 2.509571e+01, 2.458322e+01, 2.395280e+01, 2.314244e+01, 2.205922e+01, 2.051154e+01, 1.820481e+01, 1.506318e+01];
    %梯形机翼
%         Cl_dis =[2.454429e+01, 2.435579e+01, 2.405373e+01, 2.367682e+01, 2.324530e+01, 2.276954e+01, 2.225642e+01, 2.171249e+01, 2.114284e+01, 2.055110e+01, 1.993400e+01, 1.929520e+01, 1.863417e+01, 1.795095e+01, 1.725436e+01, 1.648221e+01, 1.570298e+01, 1.486106e+01, 1.390811e+01, 1.278376e+01, 1.131019e+01, 9.089243e+00];
    %11°后掠
%      Cl_dis= [2.427018e+01, 2.420620e+01, 2.415695e+01, 2.395791e+01, 2.373605e+01, 2.346473e+01, 2.318299e+01, 2.286900e+01, 2.253896e+01, 2.218677e+01, 2.182095e+01, 2.143574e+01, 2.103651e+01, 2.062482e+01, 2.020169e+01, 1.976210e+01, 1.931217e+01, 1.885467e+01, 1.838827e+01, 1.790974e+01, 1.741557e+01, 1.690650e+01, 1.637922e+01, 1.582909e+01, 1.525148e+01, 1.463911e+01, 1.398188e+01, 1.324536e+01, 1.239990e+01, 1.118425e+01, 9.442532e+00, 8.941454e+00];
     Cl_dis= [3.484079e+01, 3.466165e+01, 3.434064e+01, 3.384464e+01, 3.325592e+01, 3.259862e+01, 3.188738e+01, 3.112865e+01, 3.030735e+01, 2.944631e+01, 2.854468e+01, 2.760811e+01, 2.664938e+01, 2.565378e+01, 2.460187e+01, 2.348906e+01, 2.229025e+01, 2.097584e+01, 1.947240e+01, 1.747118e+01, 1.453580e+01, 1.228451e+01];
    Pos=linspace(0.0,1.5,22);
        
    end
    alpha0=2/180*pi;
    N=length(Pos);
    %机翼参数
    taper=0.5;
    Gama=11.3099/180*pi;
%     Gama=7.9696/180*pi;
%     Gama=0/180*pi;
    epsilon=2/3;
    % Gama=0;
    b=1.5;%半展长
    cr=0.3;%根弦长
    ct=cr*taper;
    c_ref=2/3*cr*(1+taper^2/(1+taper));%气动弦长
    S=b*cr*(1+taper)/2;
    S_cfj=b*epsilon*cr*(1+taper)/2;
    q=0.5*1.225*400;

    divide_point=(Pos(1:end-1)+Pos(2:end))/2;
    delta_b=divide_point(2:end)-divide_point(1:end-1);
    delta_b=[divide_point(1),delta_b,1.5-divide_point(end)];
    C=cr-Pos/1.5*cr*(1-taper);
    if shapes(2)==0
        Cl_dis=Cl_dis./C/(0.5*1.225*400);
    end
    
    delta_s=C.*delta_b;

    Ai=40/180*pi;
    As=50/180*pi;
    inj_loc=0.0761;
    suc_loc=0.8107;
    Bi_c=0.01;
    Bs_c=0.01*7/3;
    Bi=(Bi_c*cr+Bi_c*cr*(1-(1-taper)*epsilon))*b*epsilon/2;
    Bs=(Bs_c*cr+Bs_c*cr*(1-(1-taper)*epsilon))*b*epsilon/2;
    Angle_of_inj=atan((b*tan(Gama)+(ct-cr)*inj_loc)/b);
    Angle_of_suc=atan((b*tan(Gama)+(ct-cr)*suc_loc)/b);
    
    massrate=real_massrate*0.003/Bi;
    vinj=real_massrate/1.225/(Bi);
    Cmiu=real_massrate*vinj/(0.5*1.225*400*S_cfj);
    vsuc=real_massrate/1.225/(Bs);
    [Xl, Yl] = meshgrid(4:-2:-4, [0.03,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.135,0.15,0.16,0.17,0.18,0.2,0.22,0.25,0.275,0.3]); 
    [Xd, Yd] = meshgrid(4:-2:-4, [0.03,0.05,0.1,0.15,0.2,0.25,0.275]); 
    cl_map=table2array(readtable('cl_map_without_jet.xlsx'));
    cd_map=table2array(readtable('cd_map_without_jet.xlsx'));
    cm_map=table2array(readtable('cm_map_without_jet.xlsx'));
    cdp_map=table2array(readtable('cdp_map.xlsx'));
    cdf_map=table2array(readtable('cdf_map.xlsx'));
    Fx_map=table2array(readtable('jetFx_map.xlsx'));
    Fy_map=table2array(readtable('jetFy_map.xlsx'));
    pinj_map=table2array(readtable('pinj_map.xlsx'));
    psuc_map=table2array(readtable('psuc_map.xlsx'));
    Ptinj_map=table2array(readtable('Ptinj_map.xlsx'));
    Ptsuc_map=table2array(readtable('Ptsuc_map.xlsx'));
    alpha_map=[-6,-4,-2,0,2,4];
    Cl_map=[-0.08539475,0.104,0.296,0.485,0.667234041,0.837096513];
    Cm_map=[0.12272074,0.11815368,0.11396464,0.10980536,0.10525323,0.099945729];
    Cdp_map=[0.008149986,0.007292155,0.007278948,0.008104288,0.009828387,0.012613109];
    Cdf_map=[0.010019232,0.010385784,0.010660165,0.010845842,0.010943995,0.010955707];
    Cd_dis=zeros(N,1);
    new_Cl_dis=zeros(N,1);
    new_Cd_dis=zeros(N,1);
    Cm_dis=zeros(N,1);
    Alpha_i=zeros(N,1);
    Cd_p=zeros(N,1);
    Cd_f=zeros(N,1);
    if shapes(2)==0
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
        %         cd_eff=cdp_eff+cdf_eff;
            end
            Alpha_i(i)=alpha_i;
        end
    else
        Alpha_i=-induced_Angle/180*pi;
    end
    DeltaCl=zeros(N,1);
    DeltaCd=zeros(N,1);
    DeltaCd_p=zeros(N,1);
    DeltaCd_f=zeros(N,1);
    pinj=zeros(N,1);
    psuc=zeros(N,1);
    ptinj=zeros(N,1);
    ptsuc=zeros(N,1);
    for i=1:length(Pos)
        if Pos(i)<epsilon*1.5
     
            DeltaCl(i)=interp2(Xl, Yl, cl_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            DeltaCd_p(i)=interp2(Xd, Yd, cdp_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            DeltaCd_f(i)=interp2(Xd, Yd, cdf_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            DeltaCd(i)=DeltaCd_p(i)+DeltaCd_f(i);
            %垂直1/4弦长处
            new_Cl_dis(i)=DeltaCl(i)/cos(Alpha_i(i))-DeltaCd(i)*tan(Alpha_i(i))/cos(Alpha_i(i));
            new_Cd_dis(i)=DeltaCd(i)/cos(Alpha_i(i))+DeltaCl(i)*tan(Alpha_i(i))/cos(Alpha_i(i));
            new_Cd_dis(i)=(DeltaCd_p(i)*cos(Gama)^3+DeltaCd_f(i))/cos(Alpha_i(i))+DeltaCl(i)*tan(Alpha_i(i))/cos(Alpha_i(i))*cos(Gama)^3;
            pinj(i)=interp2(Xl, Yl, pinj_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            psuc(i)=interp2(Xl, Yl, psuc_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            ptinj(i)=interp2(Xl, Yl, Ptinj_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            ptsuc(i)=interp2(Xl, Yl, Ptsuc_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
            Cd_p(i)=DeltaCd_p(i)*cos(Gama)^3/cos(Alpha_i(i));
            Cd_f(i)=DeltaCd_f(i)/cos(Alpha_i(i));
            Cm_dis(i)=interp2(Xl, Yl, cm_map, (alpha0/cos(Gama)-Alpha_i(i))/pi*180, massrate*cos(Alpha_i(i))^2*cos(Gama-Angle_of_inj), 'linear');
        else
            Cl_eff=interp1(alpha_map,Cl_map,(alpha0/cos(Gama)-Alpha_i(i))/pi*180, 'linear');
            Cdp_eff=interp1(alpha_map,Cdp_map,(alpha0/cos(Gama)-Alpha_i(i))/pi*180, 'linear');
            Cdf_eff=interp1(alpha_map,Cdf_map,(alpha0/cos(Gama)-Alpha_i(i))/pi*180, 'linear');
            new_Cl_dis(i)=Cl_eff/cos(Alpha_i(i))-(Cdp_eff+Cdf_eff)*tan(Alpha_i(i))/cos(Alpha_i(i));
            new_Cd_dis(i)=(Cdp_eff+Cdf_eff*cos(Gama)^3)/cos(Alpha_i(i))+Cl_eff*tan(Alpha_i(i))/cos(Alpha_i(i))*cos(Gama)^3;
            Cm_dis(i)=interp1(alpha_map,Cm_map,(alpha0/cos(Gama)-Alpha_i(i))/pi*180, 'linear');
            %         Cl=Cl+Cl_dis(i);
    %         Cd=Cd+Cd_dis(i);
        end
    
    end
    Cd_p=delta_s*Cd_p/S;
    Cd_f=delta_s*Cd_f/S;
    Cl_no_jet=delta_s*new_Cl_dis/S;
    Cd_no_jet=delta_s*new_Cd_dis/S;
    bi=C*Bi_c.*delta_b;
    bs=C*Bs_c.*delta_b;
%     fprintf("Cl=%d,Cd=%d\n",Cl,Cd);
     
    DeltaP=bs*(ptinj-ptsuc)/Bs;
    fj_x=((vinj)*real_massrate)*cos(Ai-alpha0)+bi*pinj*cos(Ai-alpha0);
    fs_x=-(-(vsuc)*real_massrate)*cos(As+alpha0)-bs*psuc*cos(As+alpha0);
    fj_y=((vinj)*real_massrate)*sin(Ai-alpha0)+bi*pinj*sin(Ai-alpha0);
    fs_y=(-(vsuc)*real_massrate)*sin(As+alpha0)+bs*psuc*sin(As+alpha0);
    % Cl=Cl-interp2(Xl,Yl,Fy_map,(alpha0)/pi*180,massrate,'linear')/(q*S);
    Cl=Cl_no_jet-fj_y/(q*S)-fs_y/(q*S);
    Cd=Cd_no_jet-fj_x/(q*S)*cos(Angle_of_inj)-fs_x/(q*S)*cos(Angle_of_suc);
    Cm=delta_s*Cm_dis/S+fj_x*0.0888/(q*S)+fs_x*0.0669/(q*S)-fj_y*(0.25-0.0761)/(q*S)+fs_y*(0.25-0.8139)/(q*S);

    % Cd=Cd-interp2(Xl,Yl,Fx_map,(alpha0)/pi*180,massrate,'linear')/(0.5*1.225*400*S);
%     fprintf("Cl=%d,Cd=%d\n",Cl,Cd);
    result = [real_massrate,Cl_no_jet, Cd_no_jet, Cl, Cd,DeltaP];

    %sweep8 1.2   [5.227386e+01, 5.231627e+01, 5.201719e+01, 5.160143e+01, 5.105328e+01, 5.038221e+01, 4.960010e+01, 4.868485e+01, 4.767541e+01, 4.646466e+01, 4.503970e+01, 4.332156e+01, 4.129548e+01, 3.876133e+01, 3.400248e+01, 2.323392e+01, 2.042212e+01, 1.834275e+01, 1.655986e+01, 1.482071e+01, 1.286004e+01, 1.020453e+01];
end
 
function [Alpha_eff,Cd_p,Cd_f]=search_in_map(Cl_eff)
    Cl_map=[0.104,0.296,0.485,0.667234041,0.837096513];
    alpha_map=[-4,-2,0,2,4]/180*pi;
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
function [Alpha_eff,Cd_eff]=search_in_map_new(Cl_eff)
    Cl_map= [-0.3649565, -0.2688963, -0.1704998, -0.07037165, 0.03255308, 0.1349339, 0.2360541, 0.3365605, 0.4368884, 0.5371665, 0.6378259, 0.7388124, 0.8333122, 0.9560912, 1.04329, 1.134505, 1.225301, 1.31345, 1.396846, 1.474257, 1.544096];
    alpha_map=linspace(-10,10,21)/180*pi;
    Cd_map=[0.02028769,0.01807094,0.0165518,0.01519432,0.01452551,0.0137173,0.01319561,0.01284688,0.01282076,0.01281368,0.01282134,0.0128899,0.0125036,0.0128668,0.01367019,0.01453017,0.015516,0.01674011,0.01803182,0.01965492,0.0215665];
    N=length(Cl_map);
    for i=1:N
        if Cl_eff<Cl_map(i)
            Cd_eff=Cd_map(i-1)+(Cl_eff-Cl_map(i-1))/(Cl_map(i)-Cl_map(i-1))*(Cd_map(i)-Cd_map(i-1));
           
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

function [pos, induced_angle] = induced_angle_data_from_excel(filename)
    % 从Excel文件的第二行开始读取数据
     
    raw=readmatrix(filename);
    % 提取从第二行开始的两列数据
    data = raw(2:end, :);  % 从第二行开始
    
    % 将数据转化为数值格式并提取两列
    pos = data(:, 1);  % 第一列
    induced_angle = data(:, 2);  % 第二列
end

