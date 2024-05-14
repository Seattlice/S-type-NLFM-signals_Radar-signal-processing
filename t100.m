%S-type nonlinear FM signals - approximate solution of nonlinear FM signals using the sojourn phase principle
%S-type nonlinear FM signals_Radar signal processing
%% 清理工作区
clc;            % 清除命令行窗口的输出
clear all;      % 清除工作区的所有变量，释放内存空间

%% 非线性调频信号仿真
%% 参数设置
B = 200e6;      % 信号带宽 200 MHz
fs = 230e6;     % 采样率 230 MHz
Tp = 20e-6;     % 信号脉宽 20 μs
N = 8192;       % 点数，这里注释掉了，可能是为了适配特定的FFT算法
K = B / Tp;     % 调频斜率
ts=1/fs;        % 采样周期
t=-Tp/2:ts:Tp/2-ts;  % 时间向量，从 -Tp/2 到 Tp/2，步长 ts
f=-fs/2:fs/N:fs/2-fs/N;  % 频率向量

%% 矩形窗和时间窗函数
rectf=zeros(1,N);
for i=1:N
    if(abs(f(i))>B/2)
        rectf(i)=0;
    else
        rectf(i)=1;
    end
end   % 根据信号带宽，创建矩形窗函数

rectt=zeros(1,N);
for i=1:N
    if(i>length(t))
        rectt(i)=0;
    else
        rectt(i)=1;
    end
end   % 创建时间窗函数，确保窗函数长度与信号长度一致

%% 非线性调频信号的窗函数和群时延函数
Wf=0.54+0.46*cos(2*pi*f/B);  % Hamming 窗函数
Tf=rectf.*(Tp*f/B+0.426*Tp/pi*sin(2*pi*f/B));  % 群时延函数

%% 计算频率向量 Ft
Ft=zeros(1,N);
for k=1:length(t)
    for m=684:N
        if(t(k)<Tf(m))
            Ft(k)=(m-0.5)/N*fs-fs/2;
            break;   
        end    
    end    
end   % 通过搜索找到每个时间点对应的频率值

%% 计算累积相位和信号
fait=zeros(1,N);
for i=2:N
    fait(i)=fait(i-1)+2*pi*Ft(i)*ts;
end   % 计算累积相位

s4=exp(1i*fait);   % 计算 NLFM 信号

s4=s4.*rectt;       % 应用时间窗函数

%% 结果进行验证
T=2400;             % 验证信号的长度
s3=s4(1:4800);     % 取信号的前半部分进行分析

%% 时域和频域分析
Fc = fs;            % 载波频率，这里即为采样频率
figure;
subplot(3,1,1);
T2 = (-T+1 : T)/(Fc );  % 时域分析的时间向量
plot(T2,s3);grid on;
title('时域波形');
xlabel('时间 (秒)');
ylabel('幅度');

% 计算信号的傅里叶变换，并获取幅度谱
Y = fft(s3);        % 计算离散傅里叶变换
P2 = abs(Y/T);     % 归一化幅度谱

% 取前一半的频谱用于分析
P1 = P2(1:2*T);     
f = Fc*(0:2*(T)-1)/T;  % 构建频率向量

% 计算 dB 级别的幅度谱
Z = 20*log10(P1);

subplot(3,1,2);
plot(f/1e6, abs(fftshift(Y)),'m');grid on;
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度/db');

subplot(3,1,3);
plot(f/1e6, abs(fftshift(P1)));grid on;
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度');