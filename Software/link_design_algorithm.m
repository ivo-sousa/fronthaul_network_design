% This Link Design Algorithm is used to determine which technology 
% among Microwave Radio Trnsmission (MRT), Free Space Optics (FSO) 
% and Fiber Optics (FO) applies best to connect two points,
% namely in terms of cost.

% The MRT, FSO and FO equipment characteristics are assumed 
% to be given as inputs in the files
% MRT.dat, FSO.dat and FO.dat, respectively.

function [cheapest_cost,eq_ID] = link_design_algorithm(distance_aux,requested_debit_aux)




%Variables

% PRx [dB] -> free space propagation
% PTx [dB] -> power transmitted by the emitter
% GRx [dB] -> emitter gain
% GTx [dB] -> receiver gain
% Asis [dB] -> all losses related to equipment (OF, coaxial cables,
%receivers,...) Will be considered 3dB throughout the project.
% Ao [dB] -> Free Space attenuation
% d [km] -> Connection length
% f [GHz] -> Frequency (an array will be considered)
% lambda [m] -> Wavelength
% Aabs [dB] -> Additional atmospheric attenuation
% gamma0_o [dB/km] -> attenuation coefficient due to oxygen
% gamma0_w [dB/km] -> attenuation coefficient due to water vapour
% gamma_r [dB/km] -> attenuation due to rain
% Ri [mm/h] -> Rain intensity
% N0 [dBW] -> Thermal noise power
% b_rf [Hz] -> Effective noise bandwidth
% b0 [Hz] -> Bandwith
% B_bw -> Bandwidth difference between b0 and b_rf
% Nf [dB] -> Noise factor due to receiver





index=0;

MRT_equipment =readtable ('MRT.dat');
nr_eq_MRT=size(MRT_equipment);
FSO_equipment =readtable ('FSO.dat');
nr_eq_FSO=size(FSO_equipment);
FO_equipment =readtable ('FO.dat');
nr_eq_FO=size(FO_equipment);
 
total_nr_eq=nr_eq_MRT(1,1)+nr_eq_FSO(1,1)+nr_eq_FO(1,1);

all_eq_ref=cell(total_nr_eq,1);

scenario_data=readtable ('Scenario.dat');

switch nargin
    case 2
        d=distance_aux; %distance in km
        requested_debit=requested_debit_aux; %requested debit in Mbps
    otherwise
        d=table2array(scenario_data(1,1)); %distance in km
        requested_debit=table2array(scenario_data(1,2)); %requested debit in Mbps
end


p_time=table2array(scenario_data(1,3)); %Maximum tolerated unavailability (percentage of time)
T=table2array(scenario_data(1,4)); %Temperature
Ri_001=table2array(scenario_data(1,5)); %Rain intensity in mm/h
h_rel_aux=table2array(scenario_data(1,6)); %Relative humidity in %
h_rel=h_rel_aux/100;
ha=table2array(scenario_data(1,7)); %transmitter height in m
obs_los=table2array(scenario_data(1,8)); %height difference between tip of obstacle and line of sight in meters (negative if below)
N_fd=table2array(scenario_data(1,9)); %Number of foggy days
D_f=table2array(scenario_data(1,10));  %Duration of fog

Asis=3;




rain_frequency =readtable ('rain_k_alpha.dat');

% Earth's Atmosphere
gamma0_o=0.00613408;      % pag 58 livro STVR
gamma0_w=0.000612794;     % pag 58 livro STVR
Aabs_MRT=(gamma0_o+gamma0_w)*d;

eq_cost=ones(total_nr_eq,1);
eq_cost=eq_cost*(inf);



%Comecar o for


%%
for i=1:nr_eq_MRT(1,1) 
    index=index+1;
    
    
    debit=table2array(MRT_equipment(i,2));                    %in Mbps
    if (debit<requested_debit)
        continue % If bit rate requirement not satisfied, then skip iteration
    end
    
    all_eq_ref(index)=(table2cell(MRT_equipment(i,1)));
    f=table2array(MRT_equipment(i,3));                        %choose between (1,2,4,6,7,8,10,12,15,20,25,30GHz)
    %debit=table2array(MRT_equipment(i,3));                    %in Mbps

    PTx=table2array(MRT_equipment(i,4));                      %Transmitter Power
    GTx=table2array(MRT_equipment(i,5));
    GRx=table2array(MRT_equipment(i,6));                       
    SRx=table2array(MRT_equipment(i,7));                      %Receiver Sensitivity
    Nf=table2array(MRT_equipment(i,8));                       %Noise factor of the receiver in dB
    QAM=table2array(MRT_equipment(i,9));
    cost_fixed=table2array(MRT_equipment(i,10));
    cost_taxes=table2array(MRT_equipment(i,11));


    
    Ao=92.4+20*log10(d)+20*log10(f);

    %Rain attenuation
    
    rows_rain_frequency = rain_frequency(rain_frequency.Freq_GHz>=f,:); %loads k and alpha for assigned frequency
    k_alpha_aux=table2array(rows_rain_frequency(1,:));
    k_r=(k_alpha_aux(1,2)+k_alpha_aux(1,4))/2; %linear pol
    alpha_r=(k_alpha_aux(1,2)*k_alpha_aux(1,3) + k_alpha_aux(1,4)*k_alpha_aux(1,5))/(2*k_r);

    %if Ri_001>100
    %    d0=100;
    %else
    %    d0=35*exp(-0.015*Ri_001);
    %end
    %def=d/(1+(d/d0));
    
    gamma_r= (k_r * (Ri_001^(alpha_r))); 
    
    r_aux=1/(0.477*d^(.633)*Ri_001^(0.073*alpha_r)*f^(.123)-10.579*(1-exp(-0.024*d)));
    if r_aux>2.5
        r_aux=2.5;
    end
    def=d*r_aux;
    Ar_001=gamma_r*def;
    Ar_p= (Ar_001)* 0.12*p_time^(-0.546-0.043*log10(p_time)); %attenuation for given percentage

    %Obstacle attenuation

    r1e=17.32*sqrt(d/(4*f));       %radius of the first Fresnell zone in meters
    u=obs_los*(sqrt(2)/r1e);
    Aobs_MRT=6.9+20*log10(sqrt(((u-0.1)^2)+1)+u-0.1);

    if Aobs_MRT<0
        Aobs_MRT =0;
    end


    % signal-to-noise ratio
    PRx=PTx+GTx+GRx-Asis-Ao-Ar_p-Aobs_MRT-Aabs_MRT;          %Power detected at the receiver

    M=PRx-SRx;

    % ver pag 210

    CN_table=readtable('CN_QAM.dat');
    %QAM=1024;
    %CN = table2array(CN_table(CN_table.QAM==QAM,2));
    s = table2array(CN_table(CN_table.QAM==QAM,3));
    B_bw=0.3;
    b0=(debit*10^6)/log2(QAM);
    b_rf=(1+B_bw)*b0;
    N0=-204+10*log10(b_rf);
    CN_ipc=PRx-N0-Nf;


    sesr=0.00016;
    bber=0.000008;


    Kn=(1.4e-8)*f*d^(3.5);
    ms=8000/s;
    %Ms=10*log10(ms);

    %SESR
    mr_sesr=Kn/sesr;
    %Mr_sesr=10*log10(mr_sesr);

    mu_sesr=(mr_sesr*ms)/(ms-mr_sesr);
    Mu_sesr=10*log10(mu_sesr);

    %BBER
    alfa1=20;
    alfa2=5;
    alfa3=1;
    Nb=22500;
    rber=1e-12;
    Prber=1e-10;
    bersesr=1e-5;
    slbber=abs((log10(rber)-log10(bersesr))/(log10(Prber)-log(sesr)));

    sesr_aux=(bber-(Nb*rber)/alfa3)*(2.8*alfa2*(slbber-1)/alfa1);

    mr_bber=Kn/sesr_aux;
    %Mr_bber=10*log10(mr_bber);

    mu_bber=(mr_bber*ms)/(ms-mr_bber);
    Mu_bber=10*log10(mu_bber);


     SNRmin_sesr=CN_ipc-Mu_sesr;
     SNRmin_bber=CN_ipc-Mu_bber;

     if SNRmin_bber<SNRmin_sesr
        SNRmin=SNRmin_bber;
     else
        SNRmin=SNRmin_sesr;
     end
     
     if (SNRmin>0) && (isreal(SNRmin_bber)) && (isreal(SNRmin_sesr)) && (M>3) %&& (debit>=requested_debit)
%         X = sprintf('MRT equipment ID# %s cannot be used. SNR= %d dB & Margin= %d dB',name,SNRmin,M);
%         disp(X)

          eq_cost(index,1)=cost_fixed+cost_taxes*sqrt(d);

%         if(cost<cheapest_cost)
%             cheapest_cost=cost;
%             eq_ref=index;
%             
%         end
        
%         possible_equipment(eq_aux,1)=name;
%         possible_equipment(eq_aux,2)='MRT'
%         possible_equipment(eq_aux,3)=debit;
%         possible_equipment(eq_aux,4)=table2array(MRT_equipment(i,9));
%         possible_equipment(eq_aux,5)=0;
%         eq_aux=eq_aux+1;
        
%      else
%         X = sprintf('MRT equipment ID# %s cannot be used.',name);
%         disp(X)
          %tgt=1;
     end
 
end
 
%% Free Space Optics

p_one=(N_fd/365.25)*(D_f/24)*100;   %percentage of time (in a year) in which the visibility does not exceed 1 km
V=p_time*(1/p_one);                 %Visibility in Km

% FSO_equipment =readtable ('FSO.dat');
% nr_eq_FSO=size(FSO_equipment);


eff=1;                              
                                  % in Celsius Degrees                              

c=3e8;

h_abs=h_rel*(-0.74+90.96*exp(T/13.67)-85.4*exp(T/13.52));
w=h_abs*d*1e-3;                         %amount of water that exists in the atmosphere along distance 'd'



%Variaveis Ivo


if V<0.5
    q=0;
elseif V<1
    q=V-0.5;
elseif V<6
    q=0.16*V+0.34;
elseif V<50
    q=1.3;
else
    q=1.6;
end

%Rain Attenuation



%Arain=(0.05556+0.00848*Ri_001-3.66e-5*Ri_001^2)*10*log10(exp(1))*d;
Arain_001=d*1.076*Ri_001^.67;

Arain= (Arain_001)* 0.12*p_time^(-0.546-0.043*log10(p_time)); %attenuation for given percentage


for i=1:nr_eq_FSO(1,1)
        index=index+1;
        
        debit=table2array(FSO_equipment(i,2));              %in Mbps
        if (debit<requested_debit)
            continue % If bit rate requirement not satisfied, then skip iteration
        end
   
        all_eq_ref(index)=(table2cell(FSO_equipment(i,1)));
        lambda_nm=table2array(FSO_equipment(i,3));          %in nm
        lambda=lambda_nm*1e-9;                               %in mW
        %debit=table2array(FSO_equipment(i,3));              %in Mbps
        PTx=table2array(FSO_equipment(i,4));                %in dBW
        GTx=table2array(FSO_equipment(i,5));                %in dB
        GRx=table2array(FSO_equipment(i,6));                %in dB
        SRx=table2array(FSO_equipment(i,7));                %in dBW
        %Nf=table2array(FSO_equipment(i,8));                %Noise factor of the receiver in dB -> not necessary herein due to shot limited operation assumption
        cost=table2array(FSO_equipment(i,9));
        
   % if debit>=requested_debit
        
        opt_freq=c/lambda;
        B=debit*1e6;

        Agl=20*log10(4*pi*d)-20*log10(lambda*1e6);      %Geometric losses
        % Atmosphere Model
        if lambda == 785e-9
            Ai=0.0305;
            Ki=0.8;
            beta_i=0.112;
            Wi=54;
        else 
            Ai=0.211;
            Ki=0.802;
            beta_i=0.111;
            Wi=1.1;
        end
        if (w-Wi)>0             %alternative to using step function (w-Wi)
            flag=1;
        else
            flag=0;
        end
        
        gamma_abs=exp(-Ai*w^(1/2))*(1-flag)+Ki*((Wi/w)^beta_i)*flag;
        Aabs_FSO=gamma_abs*d;               %eq 3.9

        %Parte do Ivo

        Afog=(3.91/V)*(550/lambda_nm)^q*d;

        %No Snow
        %Turbulence

        k=(2*pi)/lambda;
        C_nsq=9.8583e-18+4.9877e-16*exp(-ha/300)+2.9228e-16*exp(-ha/1200);
        cint_var=1.23*((k)^(7/6))*C_nsq*((d*1000)^(11/6));
        Aturb=sqrt(cint_var)*2;
        
        if obs_los>0
            Aobs_FSO=+inf; %LoS obstructed -> FSO will not work
        else
            Aobs_FSO=0; %NO OBSTACLE
        end
        
        
        PRx=PTx+GTx+GRx-Asis-Agl-Aabs_FSO-Afog-Arain-Aturb-Aobs_FSO;
        M=PRx-SRx;

        SNR0=sqrt((eff*(10^((PRx+Aturb)/10)))/(B*opt_freq*2*6.63e-34));

        %SNR_av=SNR0/(sqrt((10^((Aturb)/10)))+10^(sqrt(cint_var)/10)*SNR0^2);

        %SNR_dB=10*log10(SNR0);
    %     SNR_av_dB=10*log10(SNR_av);

        BER=0.5*erfc(0.5*sqrt((SNR0)/2));

        if (BER<1e-6) && (M>3) %&& (debit>=requested_debit)
%             X = sprintf('O equipamento FSO %s pode ser usado. SNR= %d dB e Margem= %d dB',name,SNR_dB,M);
%             disp(X)
%             possible_equipment(eq_aux,1)=name;
%             possible_equipment(eq_aux,2)='FSO'
%             possible_equipment(eq_aux,3)=debit;
%             possible_equipment(eq_aux,4)=table2array(MRT_equipment(i,9));
%             possible_equipment(eq_aux,5)=0;
%             eq_aux=eq_aux+1;
%             if(cost<cheapest_cost)
%             cheapest_cost=cost;
%             eq_ref=index;
%             end
        eq_cost(index,1)=cost;
        else
%             X = sprintf('O equipamento FSO %s nao pode ser usado.',name);
%             disp(X)
        %tgt=1;
        end
%     else
%         X = sprintf('FSO equipment # %s cannot be used.',name);
%         disp(X)
%     end
   
end

%% Fiber Optics



for i=1:nr_eq_FO(1,1)
    index=index+1;
    
    BD=table2array(FO_equipment(i,3));              %in Mbps*km
    limit_debit=table2array(FO_equipment(i,2));              %in Mbps
    
    debit=min(BD/d,limit_debit);     %in Mbps
    if (debit<requested_debit)
        continue % If bit rate requirement not satisfied, then skip iteration
    end
    
    
    
    all_eq_ref(index)=(table2cell(FO_equipment(i,1)));

    Tx_min=(table2array(FO_equipment(i,4))); %Minimum transmit power dbW
    Rx_min=(table2array(FO_equipment(i,5))); %Minimum receive power dBW
    L =(table2array(FO_equipment(i,6)));     %Losses dB
    fiber_loss =(table2array(FO_equipment(i,7))); %Fiber loss dB/km
    fixed_cost=(table2array(FO_equipment(i,8)));
    costperkm=(table2array(FO_equipment(i,9)));
    
    LB=Tx_min-Rx_min;
    TL=L+d*fiber_loss;
    
    
    
    if (LB>TL+3) %&& requested_debit<=debit
        %if (cost<cheapest_cost)
        
        eq_cost(index,1)=((costperkm*d)+fixed_cost);
        
        %             cheapest_cost=cost;
        %             eq_ref=index;
        %end
    end
    
    
    
    
end
[cheapest_cost,eq_index] = min(eq_cost);
if cheapest_cost == inf
    eq_ID=[];
else
    eq_ID=all_eq_ref(eq_index);
end
end
 
