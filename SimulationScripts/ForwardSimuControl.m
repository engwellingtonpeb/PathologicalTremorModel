function [motionData] = ForwardSimuControl(SimuInfo)
%-------------------------------------------------------------------------%
%                  Federal University of Rio de Janeiro                   %
%                  Department of Biomedical Engineering                   %
%                                                                         %
%  Author: Wellington Cássio Pinheiro - MESTRADO                          %
%  Advisor: Luciano Luporini Menegaldo                                    %         
%  Date: 16/10/2018   Versão atual 07/02/2019                             %
%-------------------------------------------------------------------------%
%% Initialization
import org.opensim.modeling.*
osimModel=Model('.\models\MoBL_ARMS_tutorial_33\MoBL-ARMS OpenSim tutorial_33\ModelFiles\MoBL_ARMS_module2_4_allmuscles.osim');

%osimModel.setUseVisualizer(true);

osimState=osimModel.initSystem();

%% Model elements identification

Nstates       = osimModel.getNumStateVariables();
Ncontrols     = osimModel.getNumControls();
Ncoord        = osimModel.getNumCoordinates(); 
Nbodies       = osimModel.getNumBodies();
model_muscles = osimModel.getMuscles();
Nmuscles      = model_muscles.getSize();

% get model states
states_all = cell(Nstates,1);
for i = 1:Nstates
   states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
end

% get model muscles (controls)
Muscles = osimModel.getMuscles();  
controls_all = cell(Ncontrols,1);
for i = 1:Ncontrols
   currentMuscle = Muscles.get(i-1);
   controls_all(i,1) = cell(currentMuscle.getName());
end

% get model coordinates
Coord = osimModel.getCoordinateSet();
Coord_all = cell(Ncoord,1);
for i = 1:Ncoord
   currentCoord = Coord.get(i-1);
   Coord_all(i,1) = cell(currentCoord.getName());
end

%% Setup Joint angles 
% phini=[5,10,15,-5,-10,-15];
% psini=[85, 80, 70, 65, 60, 55];



editableCoordSet = osimModel.updCoordinateSet();
editableCoordSet.get('elv_angle').setValue(osimState, 0);
editableCoordSet.get('elv_angle').setLocked(osimState, true);

editableCoordSet.get('shoulder_elv').setValue(osimState, 0);
editableCoordSet.get('shoulder_elv').setLocked(osimState, true);

editableCoordSet.get('shoulder_rot').setValue(osimState, 0);
editableCoordSet.get('shoulder_rot').setLocked(osimState, true);

% editableCoordSet.get('elbow_flexion').setValue(osimState, deg2rad(30));
% editableCoordSet.get('elbow_flexion').setLocked(osimState, true);

%editableCoordSet.get('pro_sup').setValue(osimState, deg2rad(psini(SimuInfo.index)));
editableCoordSet.get('pro_sup').setValue(osimState, deg2rad(80));
editableCoordSet.get('pro_sup').setLocked(osimState, false);
 
editableCoordSet.get('deviation').setValue(osimState, 0);
editableCoordSet.get('deviation').setLocked(osimState, true);
 
%editableCoordSet.get('flexion').setValue(osimState, deg2rad(phini(SimuInfo.index)));
editableCoordSet.get('flexion').setValue(osimState, deg2rad(-10));
editableCoordSet.get('flexion').setLocked(osimState, false);


osimState.getY.set(40,0); %zera ativacao inicial ECRL
osimState.getY.set(42,0); %zera ativacao inicial ECRB
osimState.getY.set(44,0); %zera ativacao inicial ECU
osimState.getY.set(46,0); %zera ativacao inicial FCR
osimState.getY.set(48,0); %zera ativacao inicial FCU
osimState.getY.set(50,0); %zera ativacao inicial PQ
osimState.getY.set(52,0); %zera ativacao inicial SUP



%% Prep Simulation
%stateDerivVector = osimModel.computeStateVariableDerivatives(osimState);
osimModel.equilibrateMuscles(osimState); %solve for equilibrium similiar

%Controls function
controlsFuncHandle = @OpenSimPlantControlsFunction;
Ts=SimuInfo.Ts;
Tend=SimuInfo.Tend;
%Integrate plant using Matlab Integrator
timeSpan = [0:Ts:Tend];
integratorName = 'ode5'; %fixed step Dormand-Prince method of order 5
integratorOptions = odeset()

%% Run Simulation
set(gcf, 'color', 'white');

tic
       motionData = IntegrateOpenSimPlant(osimModel, controlsFuncHandle, timeSpan, ...
        integratorName,SimuInfo)
       
toc

%% Re-sample excitations

global U
global F
global Rq
global Dx
 %%
 a=length(U(:,1));
 timeU=0:(timeSpan(end)/a):timeSpan(end);
 timeU=timeU(1:end-1);
%%


 tsin1=timeseries(U',timeU);
 tsout=resample(tsin1,motionData.data(:,1));
 U=tsout.data;
 U=(U(:,:))';

%%
tsin2=timeseries(F((1:length(timeU)),:),timeU);
tsout=resample(tsin2,motionData.data(:,1));
F=tsout.data;

%%
tsin13=timeseries(Rq((1:length(timeU)),:),timeU);
tsout=resample(tsin13,motionData.data(:,1));
Rq=tsout.data;

%%
tsin4=timeseries(Dx((1:length(timeU)),:),timeU);
tsout=resample(tsin4,motionData.data(:,1));
Dx=tsout.data;



%% NaN
% indexNaN=find(isnan(U(:,1)));
% fillNaN=find(~isnan(U(:,1)));
% 
% for i=1:length(indexNaN)
%     
%     %F(indexNaN(i),:)=F((indexNaN(1)-1),:);
%     U(indexNaN(i),:)=U((indexNaN(1)-1),:);
%     %R(indexNaN(i),:)=R((indexNaN(1)-1),:);
% i=i+1;
% end


end

