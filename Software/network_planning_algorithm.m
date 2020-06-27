% This Network Planning Algorithm is used to optimally find 
% a fronthaul topology for wireless networks - Baseband Unit (BBU) location 
% and which technology should be used to connect to 
% each Remote Radio Head (RRH) - involving a user-defined number of points

% NOTE: it makes use of function link_design_algorithm.m - the Link Design Algorithm

function [pos_BBUs, idx_RRH_BBU, RRHs_eq_ID, network_cost] = network_planning_algorithm(nr_iterations_aux)

BBU_info=readtable('BBU.dat');

switch nargin
    case 1
        nr_iterations=nr_iterations_aux; % Replicates of K-means clustering
    otherwise
        nr_iterations=table2array(BBU_info(1,6)); % Replicates of K-means clustering
end

RRH_info=readtable('RRH.dat');
%nr_iterations=100;     %test several iterations to converge
%BBU_capacity=inf;
%BBU_cost=167000;
nr_points=size(RRH_info,1);

%BBU_capacity,BBU_cost,min_BBUs,max_BBUs
%BBU_info=readtable('BBU.dat');
RRHs_max=table2array(BBU_info(1,1)); % Max RRHs supported by a BBU
link_capacity=table2array(BBU_info(1,2)); % Max capacity of each RRH-BBU link in Mbps
BBU_cost=table2array(BBU_info(1,3)); % BBU cost in Euro
min_BBUs=table2array(BBU_info(1,4)); % minimum no. of BBUs to test
max_BBUs=table2array(BBU_info(1,5)); % maximum no. of BBUs to test

debits=table2array(RRH_info(:,3));
max_debit_link=max(debits);
if max_debit_link > link_capacity
    disp('BBUs cannot support the required bit rate of some RRH-BBU links!');
    pos_BBUs=-1;
    idx_RRH_BBU=-1;
    RRHs_eq_ID=-1;
    network_cost=inf;
    return
end

if max_BBUs > nr_points
    max_BBUs=nr_points;
end

total_sim=max_BBUs-min_BBUs+1;

pos_BBUs_sim=cell(total_sim,1);
idx_RRH_BBU_sim=zeros(nr_points,total_sim);
RRHs_eq_ID_sim=cell(nr_points,total_sim);
network_cost_sim=inf*ones(total_sim,1);


for sim_it=1:total_sim
    
    i_BBUs=sim_it+min_BBUs-1;
    
    [pos_BBUs_aux, idx_RRH_BBU_aux, RRHs_eq_ID_aux, network_cost_aux] = ...
        network_planning_algorithm_aux(RRH_info,nr_iterations,RRHs_max,BBU_cost,i_BBUs);
    if pos_BBUs_aux ~= -1 
        if network_cost_aux ~= -2
            pos_BBUs_sim(sim_it)={pos_BBUs_aux};
            idx_RRH_BBU_sim(:,sim_it)=idx_RRH_BBU_aux;
            RRHs_eq_ID_sim(:,sim_it)=RRHs_eq_ID_aux;
            network_cost_sim(sim_it)=network_cost_aux;
        end
    end
    
end

if pos_BBUs_aux == -1 % more BBUs are required
    disp('More BBUs are required!');
    pos_BBUs=-1;
    idx_RRH_BBU=-1;
    RRHs_eq_ID=-1;
    network_cost=inf;
elseif network_cost_aux == -2 % either add more BBUs or change RRHs assignment
    disp('Maximum RRHs per BBU exceeded!');
    disp('Either add more BBUs or change RRHs assignment!');
    disp('("pos_BBUs" & "idx_RRH_BBU" can be used as starting points)');
    pos_BBUs=pos_BBUs_aux;
    idx_RRH_BBU=idx_RRH_BBU_aux;
    RRHs_eq_ID=-2;
    network_cost=inf;
else
    [network_cost,best_choice]=min(network_cost_sim);
    RRHs_eq_ID=RRHs_eq_ID_sim(:,best_choice);
    idx_RRH_BBU=idx_RRH_BBU_sim(:,best_choice);
    pos_BBUs=pos_BBUs_sim{best_choice};
    
    if network_cost == inf
        disp('The equipment set does not meet the fronthaul requirements!');
        disp('Either add more BBUs or consider different equipment!');
    end
    
end

end


function [pos_BBUs, idx_RRH_BBU, RRHs_eq_ID, network_cost] = network_planning_algorithm_aux(RRH_info,nr_iterations,RRHs_max,BBU_cost,i_BBUs)

%%%%% Other Inputs
% RRH_info=readtable('RRH.dat');
% BBU_capacity=inf;
% BBU_cost=167000;

% nr_iterations=5;     %test several iterations to converge

%%%%%%

nr_points=size(RRH_info,1);

points=table2array(RRH_info(:,[1 2]));
debits=table2array(RRH_info(:,3));


if nr_points > (i_BBUs*RRHs_max) % if more RRHs than the ones supported by BBUs
    pos_BBUs=-1;
    idx_RRH_BBU=-1;
    RRHs_eq_ID=-1;
    network_cost=inf;
    
    return
end

% req_capacity=sum(debits);
% 
% k_min=ceil(req_capacity/BBU_capacity);
% if k_min > i_BBUs % not possible
%     pos_BBUs=-1;
%     idx_RRH_BBU=-1;
%     RRHs_eq_ID=-1;
%     network_cost=inf;
%     
%     return
% end



k_test=i_BBUs;

RRHs_eq_ID_temp=cell(nr_points,nr_iterations);
idx_RRH_BBU_temp=zeros(nr_points,nr_iterations);

sim_ok=0;

%while sim_ok==0

network_cost_temp=ones(nr_iterations,1);
network_cost_temp=network_cost_temp*inf;

pos_BBUs_temp=zeros(k_test,2,nr_iterations);


for j=1:nr_iterations
    [idx,C,~,D] = kmeans(points,k_test);
    
    idx_RRH_BBU_temp(:,j)=idx;
    pos_BBUs_temp(:,:,j)=C;
    D=sqrt(D)/1000;
    
    k_aux=1;
    RRHs_ok=1;
    total_cost=k_test*BBU_cost;
    
    while k_aux<=k_test && RRHs_ok==1
        idx_test= idx==k_aux;
        debits_test=debits(idx_test);
        %test_debit=sum(debits_test);
        n_test=sum(idx_test);
        %if test_debit>BBU_capacity          %associated RRHs demand more throughput than a BBU can handle
        if n_test > RRHs_max % if assigned RRHs are more than the ones supported by a BBU
            RRHs_ok=0;
            pos_BBUs_sug=C;
            idx_RRH_BBU_sug=idx;
        else
            D_test=D(idx_test,k_aux);
            %n_test=sum(idx_test);
            eq_ID_aux=cell(n_test,1);
            for dots=1:n_test
                [cost,eq_ID]=link_design_algorithm(D_test(dots),debits_test(dots));
                
                total_cost=total_cost+cost;
                if cost ~= inf
                    eq_ID_aux(dots)=eq_ID;
                end
            end
            
            RRHs_eq_ID_temp(idx_test,j)=eq_ID_aux;
            
            k_aux=k_aux+1;
        end
    end
    if RRHs_ok == 1
        sim_ok=1;
        
        network_cost_temp(j)=total_cost;

    end
    
end


%     k_test=k_test+1; %only applicable if sim_ok=0, i.e., more BBUs are required
% end

if sim_ok==1
    [network_cost,best_choice]=min(network_cost_temp);
    RRHs_eq_ID=RRHs_eq_ID_temp(:,best_choice);
    idx_RRH_BBU=idx_RRH_BBU_temp(:,best_choice);
    pos_BBUs=pos_BBUs_temp(:,:,best_choice);
else % either add more BBUs or change RRHs assignment
    pos_BBUs=pos_BBUs_sug;
    idx_RRH_BBU=idx_RRH_BBU_sug;
    RRHs_eq_ID=-2;
    network_cost=-2;
end
end

