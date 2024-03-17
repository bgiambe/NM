%-------------------------------------------------------------------------
% Code used for the normal mode decomposition analysis in the paper:
% Giambenedetti, B., Lo Bue, N., Kokoszka, F., Artale, V., and Falcini, F. 
% (2023). Multiapproach analysis of baroclinic internal tide perturbation 
% in the Ionian Sea abyssal layer (Mediterranean Sea). Geophysical Research
% Letters, 50, e2023GL104311. https://doi.org/10.1029/2023GL104311
% 
% Beatrice Giambenedetti, Last modified March 2024
%-------------------------------------------------------------------------

clear all, close all;

%---------------------------------------------------------------- Load Data

data = readtable('example_data.txt');
variables = {'Depth', 'Temperature', 'Conductivity', 'Salinity'};

if width(data)>length(variables)
    data(:,length(variables)+1:end) = [];
end
data.Properties.VariableNames = variables;

% Data Coordinates
lat= 36.3; 
long=16.1;

[Nsq2,D,Z,Nsq22, Zwork,Dwork,vect_or,d,varb,...
    varb_eig,hn,Lr,cn] = NM_fun(data, lat,...
    'SG2',[15, 1],500,1);


% ----------------------------------------------------------- Example Plots
% plot filtered and decimated N^2 profile

figure
plot(Nsq2,Z,'-k','DisplayName','Filtered Profile'), hold on
plot(Nsq22, Zwork,'.r','LineWidth',1.5, 'DisplayName','Decimated Profile')
xlabel('N^2 (s^{-2})'),ylabel('Depth (m)'),legend('location', 'best')
title('N^2')



% plot BT mode and first 3 BC

figure
    plot(vect_or(:,1),Zwork(1:end-1),'r-','LineWidth',1.5,'DisplayName','BT'), hold on
    plot(vect_or(:,2),Zwork(1:end-1),'g-', 'LineWidth',1.5,'DisplayName','BC1'), hold on
    plot(vect_or(:,3),Zwork(1:end-1),'b-','LineWidth',1.5,'DisplayName','BC2'), hold on
    plot(vect_or(:,4),Zwork(1:end-1),'c-','LineWidth',1.5,'DisplayName','BC3'), hold on
    lin = get(gca,'YLim');
    plot([0 0],lin,'--k','LineWidth',0.5)
    ylabel('Depth (m)'),xlabel('Modes')
    legend('location','best')




% plot higher modes

figure
    plot(vect_or(:,5),Zwork(1:end-1),'r-','LineWidth',1.5,'DisplayName','BC4'), hold on
    plot(vect_or(:,6),Zwork(1:end-1),'g-', 'LineWidth',1.5,'DisplayName','BC5'), hold on
    plot(vect_or(:,7),Zwork(1:end-1),'b-','LineWidth',1.5,'DisplayName','BC6'),hold on
    plot(vect_or(:,8),Zwork(1:end-1),'c-','LineWidth',1.5,'DisplayName','BC7')
    lin = get(gca,'YLim');
    plot([0 0],lin,'--k','LineWidth',0.5)
    ylabel('Depth (m)'),xlabel('Modes')








