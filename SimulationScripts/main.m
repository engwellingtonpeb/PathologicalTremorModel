%-------------------------------------------------------------------------%
%                  Federal University of Rio de Janeiro                   %
%                  Department of Biomedical Engineering                   %
%                                                                         %
%  Author: Wellington Cássio Pinheiro - MsC Requirement                   %
%  Advisor: Luciano Luporini Menegaldo                                    %         
%  Date:01/11/20                                                      %
%-------------------------------------------------------------------------%
%                                                                         %
% 1) This code run a OpenSim forward dynamics simulation to produce a     %
% wrist motion dataset.                                                   %
%                                                                         %
% 2) A linear state space model is estimatimated forW1=w1 activation,          %
% contraction and joint dynamics.                                         %
%                                                                         %
% 3) Model validation routines are performed                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

%close all
run=true;

SimuInfo=struct %information about simulation parameters
MotionValida=struct;
MotionValida.data=[];
SimuInfo.estimation=false;
SimuInfo.validation=false;
SimuInfo.ControlTest=false;
while (run)
%% Variable Definition (DIALOG)
    disp('===============================================================')
    disp('(1) - Estimar Modelo')
    disp('(2)-  Validar Modelo')
    disp('(3)-  Control Test')
    opt1=input('Digite a opção:','s')

    if (opt1=='1')
     SimuInfo.estimation=true;
     SimuInfo.Tend=input('Digite o tempo simulado em segundos:')
    end
    if(opt1=='2')
     SimuInfo.validation=true; 
    end
    
    if(opt1=='3')
     SimuInfo.ControlTest=true; 
     SimuInfo.Tend=input('Digite o tempo simulado em segundos:')
    end

    if (opt1=='2');
      disp('Qual teste deseja realizar?')
      disp('(1) - Variar Amplitude')
      disp('(2) - Variar Frequência')
      disp('(3) - Variar Ruído')
      disp('(4) - Variar Seq. de recrutamento')
      disp('(5) - Modelo em Cascata')
      opt2=input('Digite a opção:','s')
      opt=strcat(opt1,opt2);
    else
      opt=opt1;
    end

    SimuInfo.opt=opt;
    %% Generate Identification Dataset
while (SimuInfo.estimation)
   

    SimuInfo.Amplitude=0.20;
    SimuInfo.Freq=50;
    %SimuInfo.Noise=40;
   
    tic
       [motionData,U,F,R]=ForwardGetMotion(SimuInfo); % Generate model train dataset
    toc 

%%    
ShowMotionData(motionData,U,F)   

%% Get Torque 

TAU_f=F(:,(1:7)).*R(:,(1:7)); % Torque of Flexion (positive) / Extension (negative) of each muscle
TAU_WristFlex=sum(TAU_f,2);

TAU_p=F(:,(1:7)).*R(:,(8:14)); % Torque of Pronation (positive) / Suppination (negative) of each muscle
TAU_WristPro=sum(TAU_p,2);



%% saving simulation

string1=num2str(SimuInfo.Amplitude);
string1=strrep(string1,'.','');
string1=strcat('A',string1);
string1=strcat(string1,'_');
string2=num2str(SimuInfo.Freq);
string2=strcat('F',string2);
string2=strcat(string2,'_');
string1=strcat(string1,string2);
date = datestr(datetime('now')); 
date=regexprep(date, '\s', '_');
date=strrep(date,':','_');
date=strrep(date,'-','_');
date=strcat(date,'_');
string1=strcat(date,string1);
string1=strcat(string1,'_Motion')

save(string1,'motionData','U','F','R','TAU_f','TAU_p')

%% Call Idtf function
% FITa01=[];
% FITc01=[];
% 
% FITac02=[];
% 
% 
% FITac03=[];
% FITac04=[];

FIT=[];


%FITj=[];

tic
for i=1:50    %Loop to verify order with better fit
    nx1=i+1;
    nx2=i+1;
    nx3=i+1;
    disp(nx1)

    [fit] = MimoIdtf(motionData,U,TAU_p,nx1,nx2,nx3)

    FIT=[FIT,fit];
   

end
toc
%%
disp('Melhor ordem de modelos')
%ativação
    [minNumCol,minIndexLin]=min(FITa); %Among minimum fit for multiple channels and order, find the maximum fit
    [minNum,cola]=max(minNumCol);
    minNum
    cola+1
%contração
    [minNumCol,minIndexLin]=min(FITc);
    [minNum,colc]=max(minNumCol);
     minNum
     colc+1
%posição
    [minNumCol,minIndexLin]=min(FITj);
    [minNum,colj]=max(minNumCol);
    minNum
    colj+1
%%
%     nx1=cola+1; %ajuste pq começa com sistemas de ordem 2
%     nx2=colc+1;
%     nx3=colj+1;
    [Sys_ativ,ativ,Sys_cont,cont,Sys_joint,joint,fita,fitc,fitj] = MimoIdtf(motionData,U,F,nx1,nx2,nx3);


    %% Individual block check
    ativ1=ativ;
    cont1=cont;
    joint1=joint;
    figure
    compare(ativ1,Sys_ativ);
    figure
    compare(cont1,Sys_cont);
    figure
    compare(joint1,Sys_joint);
    %% Cascade system check
uecrl=U(1,:);
% uecrb=U(2,:);
% uecu=U(3,:);
% ufcr=U(4,:);
ufcu=U(5,:);
%upq=zeros(1,length(ufcu));
%     ativ=iddata([motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],[uecrl',uecrb',uecu',ufcr',ufcu'],0.001);
%     ativ.InputName={'uECRL';'uECRB';'uECU';'uFCR';'uFCU'};
%     ativ.OutputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};

    ativ=iddata([motionData.data(:,42),motionData.data(:,50)],[uecrl',ufcu'],0.0001);
    ativ.InputName={'uECRL';'uFCU'};
    ativ.OutputName={'ECRL Activation';'FCU Activation'};
    [ya,fita,x0a]=compare(ativ,Sys_ativ);
    figure
    compare(ativ,Sys_ativ);
    ya=ya.OutputData;
%%
%     cont=iddata([motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],[ya(:,1),ya(:,2),ya(:,3),ya(:,4),ya(:,5)],0.001);
%     cont.InputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};
%     cont.OutputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length'};
   
    F1=F(:,1);
    F5=F(:,5);
    cont=iddata([-F1,F5],[ya(:,1),ya(:,2)],0.0001);
    cont.InputName={'ECRL Activation';'FCU Activation'};
    cont.OutputName={'ECRL Force';'FCU Force'};
    [yc,fitc,x0c]=compare(cont,Sys_cont);
    figure
    compare(cont,Sys_cont);
    yc=yc.OutputData;
    yc=yc;
%%

%     joint=iddata([rad2deg(motionData.data(:,19)),rad2deg(motionData.data(:,39)),rad2deg(motionData.data(:,17)),rad2deg(motionData.data(:,37))],[yc(:,1),yc(:,2),yc(:,3),yc(:,4),yc(:,5)],0.001);
%     joint.InputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length';};
%     joint.OutputName={'theta';'dtheta';'phi';'dphi'}; 
    figure
    
    joint=iddata([(motionData.data(:,19)),(motionData.data(:,17))],[-yc(:,1), yc(:,2)],0.0001);
    joint.InputName={'ECRL Force';'FCU Force'};
    joint.OutputName={'theta';'phi'};

    compare(joint,Sys_joint);
    [yj,fitj,x0j]=compare(joint,Sys_joint,5);
   
%%    
  
estimation=false;
break;
end

    %% Validation
while(SimuInfo.validation)
    
    disp('(1) - Usar modelo que está no workspace')
    disp('(2) - Carregar modelo')
    
    %% Call Validation Test#01 
    tic
    if(opt2=='1')
    SimuInfo.Tend=input('Digite o tempo simulado em segundos:')
    SimuInfo.Amplitude=0.05
    SimuInfo.Freq=10;
    SimuInfo.Noise=400;
    load('Modelos_04_07.mat')
    Fita=[];
    Fitc=[];
    Fitj=[];
    for i=1:20

        
        [motionData,U]=ForwardGetMotionStiff(SimuInfo); 
        i=i+1;
        
        ativ=iddata([motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],[U(1,:)',U(2,:)',U(3,:)',U(4,:)',U(5,:)'],0.001);
        ativ.InputName={'uECRL';'uECRB';'uECU';'uFCR';'uFCU'};
        ativ.OutputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};

        
        cont=iddata([motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],[motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],0.001);
        cont.InputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};
        cont.OutputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length'};
     
        
        joint=iddata([motionData.data(:,19),motionData.data(:,39),motionData.data(:,17),motionData.data(:,37)],[motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],0.001);
        joint.InputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length'};
        joint.OutputName={'theta';'dtheta';'phi';'dphi'};
        
        MotionValida.data=[MotionValida.data;motionData.data];
        SimuInfo.Amplitude=(SimuInfo.Amplitude)+0.05;
       
        figure
        [y,fita,x]=compare(ativ,Sys_ativ);
        figure
        [y,fitc,x]=compare(cont,Sys_cont);
        figure
        [y,fitj,x]=compare(joint,Sys_joint);
        
        FITa=[FITa; fita];
        FITc=[FITc; fitc];
        FITj=[FITj; fitj];
        
    end
    
    toc
        break;

    end
    
   
    %% Call Validation Test#02
tic
    if(opt2=='2')
    SimuInfo.Tend=input('Digite o tempo simulado em segundos:')
    SimuInfo.Amplitude=0.05
    SimuInfo.Freq=2;
    SimuInfo.Noise=40;
    FITa=[];
    FITc=[];
    FITj=[];
    load('Modelos_04_07.mat')
    for i=1:20
 
        [motionData,U]=ForwardGetMotion(SimuInfo); 
        i=i+1;
        
        ativ=iddata([motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],[U(1,:)',U(2,:)',U(3,:)',U(4,:)',U(5,:)'],0.001);
        ativ.InputName={'uECRL';'uECRB';'uECU';'uFCR';'uFCU'};
        ativ.OutputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};

        
        cont=iddata([motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],[motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],0.001);
        cont.InputName={'ECRL Activation';'ECRB Activation';'ECU Activation';'FCR Activation';'FCU Activation'};
        cont.OutputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length'};
     
        
        joint=iddata([motionData.data(:,19),motionData.data(:,39),motionData.data(:,17),motionData.data(:,37)],[motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],0.001);
        joint.InputName={'ECRL Fiber Length';'ECRB Fiber Length';'ECU Fiber Length';'FCR Fiber Length';'FCU Fiber Length'};
        joint.OutputName={'theta';'dtheta';'phi';'dphi'};
        
        MotionValida.data=[MotionValida.data;motionData.data];
        SimuInfo.Freq=(SimuInfo.Freq)+2;
       
        figure
        [ya,fita,x0a]=compare(ativ,Sys_ativ);
        compare(ativ,Sys_ativ);
        figure
        [yc,fitc,x0c]=compare(cont,Sys_cont);
        compare(cont,Sys_cont);
        figure
        [yj,fitj,x0j]=compare(joint,Sys_joint);
        compare(joint,Sys_joint);
        
        

        
        FITa=[FITa, fita];
        FITc=[FITc, fitc];
        FITj=[FITj, fitj];
        
    toc    
    end
    
    toc
    break;
    end
    %% Call Validation Test#03
    if(opt2=='3')
    SimuInfo.Tend=input('Digite o tempo simulado em segundos:')
    SimuInfo.Amplitude=0.1;
    SimuInfo.Freq=2;
    SimuInfo.Noise=100;
    load('Modelos_04_07.mat')
    for i=1:10
        [motionData,U]=ForwardGetMotion(SimuInfo); 
        ativ=iddata([motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],[U(1,:)',U(2,:)',U(3,:)',U(4,:)',U(5,:)'],0.001);
        cont=iddata([motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],[motionData.data(:,42),motionData.data(:,44),motionData.data(:,46),motionData.data(:,48),motionData.data(:,50)],0.001);
        joint=iddata([motionData.data(:,19),motionData.data(:,39),motionData.data(:,17),motionData.data(:,37)],[motionData.data(:,43),motionData.data(:,45),motionData.data(:,47),motionData.data(:,49),motionData.data(:,51)],0.001);
        MotionValida.data=[MotionValida.data;motionData.data];
        SimuInfo.Noise=(SimuInfo.Noise)-6;
        
        figure
        compare(ativ,Sys_ativ);
        figure
        compare(cont,Sys_cont);
        figure
        compare(joint,Sys_joint);
    end
    end
    break;

    %% Call Validation Test#04
    if(opt2=='4')

    end
    break;
end
    
while(SimuInfo.ControlTest)

    for k=1:1
        clearvars -except SimuInfo k
  
        load('26_Nov_2018_08_43_40_A0075_F50_ControllerDiscrete.mat')

        %Distribuição de um paciente especíico
        %load('distrib_tremor_paciente01.mat') % paciente
        
        %Distribuição Genérica de frequências do tremor
        X = makedist('Normal','mu',5.6,'sigma',1);%generico
        P=[];
        N=570;
        for i=1:N
          P(i)=random(X,1,1);
        end
        
%         SimuInfo.Kz=Kz;
        SimuInfo.Ts=.0001;
        SimuInfo.Kz=c2d(K,SimuInfo.Ts)
        
        [Ak,Bk,Ck,Dk]=ssdata(SimuInfo.Kz);
        
        SimuInfo.Ak=Ak;
        SimuInfo.Bk=Bk;
        SimuInfo.Ck=Ck;
        SimuInfo.Dk=Dk;
        
        SimuInfo.P=P;
        pd = makedist('Uniform','lower',1,'upper',length(P));
        SimuInfo.pd=pd;
        

        
        PhiRef=0;%makedist('Normal','mu',0,'sigma',4);
        PsiRef=80;%makedist('Normal','mu',60,'sigma',0);
        
        SimuInfo.PhiRef=PhiRef;
        SimuInfo.PsiRef=PsiRef;
        SimuInfo.index=k;


        tic
        [motionData]=ForwardSimuControl(SimuInfo); % Generate model train dataset
        
        Data = load('handel.mat');
        sound(Data.y, Data.Fs)
        toc 
        run=false;


        formatOut = 'yyyy/mm/dd/HH/MM/SS';
        date=datestr(now,formatOut);
        date=strrep(date,'/','_');

        indir=pwd;
        indir=strcat(indir,'\Simulações\TremorPadrao');
        filename=strcat(date,'_MScPaperKL_COE731')
        extension='.mat';
        motionFilename=fullfile(indir,[filename extension]);

        global U
        global F
        global Rq
        global Dx
        save(motionFilename,'motionData','SimuInfo','U','F','Rq','Dx');
    k=k+1;
    end

    break;
end
end

     Data = load('handel.mat');
     sound(Data.y, Data.Fs)