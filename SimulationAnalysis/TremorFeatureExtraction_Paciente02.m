%-------------------------------------------------------------------------%
%                  Federal University of Rio de Janeiro                   %
%                  Department of Biomedical Engineering                   %
%                                                                         %
%  Author: Wellington Cássio Pinheiro, MSc.                               %
%  Advisor: Luciano Luporini Menegaldo                                    %         
%  Date: 15/10/2020                                                       %
%-------------------------------------------------------------------------%
%
% Analysis of tremor signals from voluntary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

load('p02_59_f_et_LA_20yr.mat')

%% load signal
add=1;
P=[];
P1=[];
%while (add)
    
      

     pd011=table2array(p02_2); % no artigo foram utilizados os dados etmaria2
    

     %EMG Extensores e Flexores
    Fs_emg=1925.93;
    t_emg=(pd011(:,1));
    N=length(t_emg);% comprimento do vetor de tempo
    ts=1/Fs_emg; % intervalo de tempo entre duas amostras
    Ttotal=ts*(N-1);
    t_emg=[0:ts:Ttotal];
    Xemg=([pd011(:,2) pd011(:,10)]); %Xemg(:,1)= Flexor  || Xemg(:,2)= Extensor

    %ACC Dorso da Mão
    Fs_acc=148.148;
    ts=1/Fs_acc; % intervalo de tempo entre duas amostras
    t_acc=[0:ts:Ttotal];
    N=length(t_acc);
    Xacc=([pd011((1:N),20) pd011((1:N),22) pd011((1:N),24)]); %Xacc=[Acc_X Acc_Y Acc_Z]
    % 
    %Gyro Dorso da Mão
    Fs_gyro=148.148;
    % comprimento do vetor de tempo
    ts=1/Fs_gyro; % intervalo de tempo entre duas amostras

    % Ttotal=ts*(N-1);
    t_gyro=[0:ts:Ttotal];
    N=length(t_gyro);
    Xgyro=([pd011((1:N),26) pd011((1:N),28) pd011((1:N),30)]); %Xacc=[Gyro_X(pro_sup) Gyro_Y(flex) Gyro_Z(rad_ulnar)]

     
     
     
     
    ts_gyro=1/Fs_gyro;
    Tjan=1;
    Njan=round(Tjan/ts_gyro); %qtd ptos na janela
    r=rectwin(Njan);%Define janela RETANGULAR de comprimento Njan
    h=hamming(Njan);%Define janela HAMMING de comprimento Njan
    N=length(t_gyro);
    w1=(floor(N/6));
    
    
        fc=15;
    delta=10;
    Wp=(fc-delta)/(Fs_gyro/2);
    Ws=(fc+delta)/(Fs_gyro/2);
    Rp=0.1;
    Rs=60;
    [Ng,Wn] = buttord(Wp, Ws, Rp, Rs);
    [B,A] = butter(Ng,Wn);
    
    Xgyro_f(:,1)=filtfilt(B,A,Xgyro(:,1));
    Xgyro_f(:,2)=filtfilt(B,A,Xgyro(:,2));
    Xgyro_f(:,3)=filtfilt(B,A,Xgyro(:,3));


    
 for ij=1:6 %6 JANELAS DE 10 SEGUNDOS
    
    F=0:.1:20;
    overlap=.5*Njan; % 50% overlap
    [s,w,t] =spectrogram(Xgyro_f(((w1*ij-w1+1):(w1*(ij+1)-w1-1)),2),h,overlap,F,Fs_gyro,'yaxis');
    s=abs((s)); %(ANALISE DE JANELAS DE 10 SEGUNDOS)
    s=s./max(max(s)); %normaliza a amplitude (q nao é importante na analise)
    figure(ij)
    surf( t, w, s );
%     title('Espectrograma s/ Overlap - Janela Hamming')
    ylabel('Frequência(Hz)')
    xlabel('Tempo(s)')
    zlabel('Amplitude')
    colormap jet

     
     
     %FREQUENCY HISTOGRAM

        [k,l]=size(s);
        
        for i=1:l
            [val,k]=max(s(:,i));
            P=[P F(k)];
        end
  %% Limit Cycle
        Xbase(:,1)=Xgyro(:,1)-mean(Xgyro(:,1));
        Xbase(:,2)=Xgyro(:,2)-mean(Xgyro(:,2));
        Xbase(:,3)=Xgyro(:,3)-mean(Xgyro(:,3));
        
X=[cumtrapz(t_gyro,Xbase(:,1)) cumtrapz(t_gyro,Xbase(:,2)) cumtrapz(t_gyro,Xbase(:,3))];
        
  
% 
%  %% EMG X Ativação
%  
%  %Low Pass Filtering
 Xemg(:,1)=Xemg(:,1)./max(Xemg(:,1));
 Xemg(:,2)=Xemg(:,2)./max(Xemg(:,2));
%  load('filter.mat')
 
 
 end    
     
 
Phi_ref=X(:,2);
Phidot_ref=Xgyro_f(:,2);
Psi_ref=X(:,1);
Psidot_ref=Xgyro_f(:,1);

fc=15;
delta=10;
Wp=(fc-delta)/(Fs_emg/2);
Ws=(fc+delta)/(Fs_emg/2);
Rp=0.1;
Rs=60;
[Ng,Wn] = buttord(Wp, Ws, Rp, Rs);
[B,A] = butter(Ng,Wn);
%%
y=filter(B,A,Xemg(:,1));



%% Sinais do Modelo

t_simu=[];
Phi_simu=[];
Psi_simu=[];
Phidot_simu=[];
Psidot_simu=[];
a_ecrl=[];
a_fcu=[];

% load('2020_12_13_21_14_07_MScPaperKL_paciente02.mat')
load('2020_12_10_11_06_01_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

% load('2020_12_13_21_14_07_MScPaperKL_paciente02.mat')
load('2020_12_10_11_55_55_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

% load('2020_12_13_22_30_13_MScPaperKL_paciente02.mat')
load('2020_12_10_12_45_49_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

% load('2020_12_13_23_46_21_MScPaperKL_paciente02.mat')
load('2020_12_10_14_14_57_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

% load('2020_12_14_01_02_34_MScPaperKL_paciente02.mat')
load('2020_12_10_15_05_07_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

% load('2020_12_14_02_18_53_MScPaperKL_paciente02.mat')
load('2020_12_10_15_54_29_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU
% 
load('2020_12_10_16_43_59_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

load('2020_12_10_17_33_37_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((10000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((10000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((10000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((10000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

load('2020_12_10_18_23_08_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((40000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((40000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((40000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((40000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

load('2020_12_14_07_18_57_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((40000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((40000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((40000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((40000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

load('2020_12_14_08_33_12_MScPaperKL_paciente02.mat')

Phi_simu=[Phi_simu;rad2deg(motionData.data((40000:end),19))]; %~10000 é o numero da amostra onde entra o oscilador
Psi_simu=[Psi_simu;rad2deg(motionData.data((40000:end),17))];
Phidot_simu=[Phidot_simu;rad2deg(motionData.data((40000:end),39))];
Psidot_simu=[Psidot_simu;rad2deg(motionData.data((40000:end),37))];
a_ecrl=[a_ecrl; motionData.data((10000:end),44)];  %ativ ECRL
a_fcu=[a_fcu; motionData.data((10000:end),52)];  %ativ FCU

Ts=.0001;
t_simu=0:Ts:(length(Phi_simu)-1)*Ts;

[Phi_simu,t_new] = resample(Phi_simu,t_simu,Fs_gyro);
[Psi_simu,t_new] = resample(Psi_simu,t_simu,Fs_gyro);

Psi_simu=Psi_simu-mean(Psi_simu);

[Phidot_simu,t_new] = resample(Phidot_simu,t_simu,Fs_gyro);
[Psidot_simu,t_new] = resample(Psidot_simu,t_simu,Fs_gyro);

Phi_simu=Phi_simu(825:end);
Psi_simu=Psi_simu(825:end);
Phidot_simu=Phidot_simu(825:end);
Psidot_simu=Psidot_simu(825:end);

[a_ecrl,t_emg_simu] = resample(a_ecrl(1: length(t_simu)),t_simu,Fs_emg);
[a_fcu,t_emg_simu] = resample(a_fcu(1: length(t_simu)),t_simu,Fs_emg);

%spectrogram from simulation


for ij=1:6 %6 JANELAS DE 10 SEGUNDOS
    
    F=0:.1:20;
    overlap=.5*Njan; % 50% overlap
    [s,w,t] =spectrogram(Phidot_simu(((w1*ij-w1+1):(w1*(ij+1)-w1-1)),1),h,overlap,F,Fs_gyro,'yaxis');
    s=abs((s)); %(ANALISE DE JANELAS DE 10 SEGUNDOS)
    s=s./max(max(s)); %normaliza a amplitude (q nao é importante na analise)
    figure(6+ij)
    surf( t, w, s );
%     title('Espectrograma s/ Overlap - Janela Hamming')
    ylabel('Frequência(Hz)')
    xlabel('Tempo(s)')
    zlabel('Amplitude')
    colormap jet

     
     
     %FREQUENCY HISTOGRAM

        [k,l]=size(s);
        
        for i=1:l
            [val,k]=max(s(:,i));
            P1=[P1 F(k)];
        end

end       
figure
histogram(P)
hold on
histogram(P1,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','-.','LineWidth',1.5)

%%
[Table] = MetricsTable(P,P1,Phi_ref,Phi_simu,Psi_ref,Psi_simu,Phidot_ref,Phidot_simu,Psidot_ref,Psidot_simu)

%% Ciclos

figure

plot(Phi_ref,Phidot_ref,'b','LineWidth',2)
grid 
hold on
plot(Phi_simu,Phidot_simu,'r-.','LineWidth',.5)

figure
plot(Psi_ref,Psidot_ref,'b','LineWidth',2)
grid 
hold on
plot(Psi_simu,Psidot_simu,'r-.','LineWidth',.5)


figure
histogram(Phi_ref)
hold on
histogram(Phi_simu,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','-.','LineWidth',1.5)


figure
histogram(Psi_ref)
hold on
histogram(Psi_simu,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','-.','LineWidth',1.5)


figure
histogram(Phidot_ref)
hold on
histogram(Phidot_simu,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','-.','LineWidth',1.5)


figure
histogram(Psidot_ref)
hold on
histogram(Psidot_simu,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','-.','LineWidth',1.5)
