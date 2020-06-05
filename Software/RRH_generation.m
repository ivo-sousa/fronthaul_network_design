% This script generates random RRH positions and associated debits 
% (to be used to extract conclusions by using network_planning_algorithm.m - the Network Planning Algorithm)

%save aux_workspace
%clearvars

l1=20;         %                            % Test area length (km)
l2=20;         %                            % Test area width (km)
density=0.015;  %                            % RRH density by km^2


nr_points=ceil(density*(l1*l2));       % RRH positions to be generated

X=floor(rand(nr_points,1)*1000*l1);
Y=floor(rand(nr_points,1)*1000*l2);

possible_debits=1000;
%possible_debits=[100,400,600,1000,1200,3000,6000,2000]; % Debits to be considerated (random choice among these)

points=[X,Y];

debits=zeros(nr_points,1);
for i=1:nr_points
    debits(i)=possible_debits(ceil(rand*size(possible_debits,2)));
end

x=[l1*500,l2*500];

% figure
% hold on
% plot(points(:,1),points(:,2),'*');
% plot(x(1,1),x(1,2),'x');

%save RRH.mat

T = table(X,Y,debits,'VariableNames',{'Position_X_m' 'Position_y_m' 'Debit_Mbps'});
writetable(T,'RRH.dat')

%load aux_workspace