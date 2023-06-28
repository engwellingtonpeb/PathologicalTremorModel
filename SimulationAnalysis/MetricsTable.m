function [Table] = MetricsTable(P,P1,Phi_ref,Phi_simu,Psi_ref,Psi_simu,Phidot_ref,Phidot_simu,Psidot_ref,Psidot_simu)
close all
clc
%% freq
w=2*iqr(P1)*length(P1)^(-1/3);
edges=[min(P1),max(P1)];

figure(1)
[Metrics] = ModelMetrics(P,P1,edges,w);

col1=[Metrics.dKL;Metrics.dI; abs(Metrics.centroid_P(1)-Metrics.centroid_Q(1)); Metrics.JSD]


%% phi
w=2*iqr(Phi_ref)*length(Phi_ref)^(-1/3);
edges=[min(Phi_simu),max(Phi_simu)];

figure(2)
[Metrics] = ModelMetrics(Phi_ref,Phi_simu,edges,w);

col2=[Metrics.dKL;Metrics.dI; abs(Metrics.centroid_P(1)-Metrics.centroid_Q(1)); Metrics.JSD]

%% psi
w=2*iqr(Psi_simu)*length(Psi_simu)^(-1/3);
edges=[min(Psi_simu),max(Psi_simu)];

figure(3)

[Metrics] = ModelMetrics(Psi_ref,Psi_simu,edges,w);

col3=[Metrics.dKL;Metrics.dI; abs(Metrics.centroid_P(1)-Metrics.centroid_Q(1)); Metrics.JSD]
%% phidot
w=2*iqr(Phidot_simu)*length(Phidot_simu)^(-1/3);
edges=[min(Phidot_simu),max(Phidot_simu)];

figure(4)
[Metrics] = ModelMetrics(Phidot_ref,Phidot_simu,edges,w);
col4=[Metrics.dKL;Metrics.dI; abs(Metrics.centroid_P(1)-Metrics.centroid_Q(1)); Metrics.JSD]
%% psidot
w=2*iqr(Psidot_simu)*length(Psidot_simu)^(-1/3);
edges=[min(Psidot_simu),max(Psidot_simu)];

figure(5)
[Metrics] = ModelMetrics(Psidot_ref,Psidot_simu,edges,w);
col5=[Metrics.dKL;Metrics.dI; abs(Metrics.centroid_P(1)-Metrics.centroid_Q(1)); Metrics.JSD]


Table=[col1 col2 col3 col4 col5]
end

