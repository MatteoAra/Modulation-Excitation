%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THE FOLLOWING LINE IN CASE OF DIFFERENT INPUT FILE...
clear all; close all 
FILENAME='34570_Pd_Al2O3_H2xO2_300C_crop_rebin_crop_FFT_k12_cut_2879.dat';     %Remember to remove the commented header in the .dat file!

phistart=0;
phistop=350;
phistep=10;

K=1;            % User-defined K value set to 1
timestep=0.5;   % Time step (s) between acquisition of different spectra during the experiment
numframes=144;  % Number of frames in single cycle (144 for H2xO2 and 288 for CH4_CO2 CO)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%import_data=load(FILENAME);
%import_data=aa(:,2:end);
data_buffer=load(FILENAME); 
energies=data_buffer(:,1);
%import_data=period_merger(data_buffer(:,2:end), numframes); % summing 16 cycles 
import_data=period_merger(data_buffer(:,2:145), numframes); % 1st cycle
%only, change if H2xO2 (145)

numscans=size(import_data,2);    % Number of scans over the desired integration period T
phival=deg2rad([phistart:phistep:phistop]);
timestamp=[1:numscans].*timestep;
T=numscans.*timestep;

datatimesin=import_data.*sin((K*2*pi./T).*timestamp);
datatimecos=import_data.*cos((K*2*pi./T).*timestamp);

intI=[];
for index2=1:length(phival)
    for index1=1:length(import_data)
        buf=(datatimesin(index1,:)*cos(phival(index2)))+(datatimecos(index1,:)*sin(phival(index2)));
        %buf = import_data(index1,:) .* sin( (K*2*pi./T).*timestamp + phival(index2) );
        intI(index1,index2)=(2/T).*trapz(buf);
    end
end


if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical post-processing
plot(energies,intI,'Linewidth',1.5); 
set(gca,'Fontsize',18); xlabel(['R (1/' char(197),')']); ylabel('Intensity (Arb. Units)'); 

disp('Now click where you would like to view the demodulation analysis and then press return')
[Xsel,Ysel] = ginput;
X_index=[];
% Definition of legend
Legend=cell(length(Xsel),1);
for ind=1:length(Xsel)
    Legend{ind}=string(['R = ', num2str(Xsel(ind))]);
end
% Plot generation
figure(2)
B=ones(size(phival,2),length(Xsel)); 
for ind = 1:length(Xsel)
    [c X_index(ind)] = min(abs(energies-Xsel(ind)));
    B(:,ind)=intI(X_index(ind),:);
    subplot(2,1,1)
    plot(phival*(180/pi),intI(X_index(ind),:)); hold on
end
B=[((180/pi).*phival)' B];
TITLE=['Demodulation angle for the ', num2str(length(Xsel)),' user-defined R value(s)']; 
title(TITLE);
set(gca,'Fontsize',18,'box','on'); xlabel(['Demodulation Angle (\circ)']); ylabel('Intensity (Arb. Units)'); xlim([0 350]);
ll=legend(Legend); set(ll,'box','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold analysis
%thrs=input('Define Signal Threshold (from 0 to 1, corrisp to 0% to 100%): ')
thrs=0.33; % Only allowing values from 0 to 1 (namely 0% to 100%)

if thrs<0
    thrs=0;
elseif thrs>1
    thrs=1;
end

C=[];
refI=thrs*max(max(abs(intI)));
for ind=1:size(intI,2)
    indexAbove=find(abs(intI(:,ind)) > refI);     
    subplot(2,1,2)
    if length(indexAbove)>0
        temparray=ones(length(indexAbove),1).*(phival(ind)*(180/pi));
        buff=[energies(indexAbove),temparray]; C=[C; buff];
        plot(energies(indexAbove),temparray,'bo' ); hold on
    end
end

TITLE=['Phase Value above Absolute Intensity threshold of ', num2str(thrs*100),'%']; 
title(TITLE);
set(gca,'Fontsize',18,'box','on'); xlabel(['R (1/' char(197),')']); ylabel(['Phase Value (\circ)']); 
xlim([min(energies) max(energies)]); ylim([-5 365]);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/O 

ioval=input('Should I save the data? (y/n)', 's');
switch ioval
    case 'y'
        ifval=1;
    case 'n'
        ifval=0;
end
        
if ifval  % change to 1 to save in ascii format demodulation-processed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/O - Saving parameters here below 
OUT_FILENAME=[FILENAME,'_PhaseDemod.dat'];
A=[energies intI];
eval(string(strcat('save',{' '}, OUT_FILENAME,' A -ascii')));

OUT_FILENAME=[FILENAME,'_IvsPhase.dat'];
eval(string(strcat('save',{' '}, OUT_FILENAME,' B -ascii')));    

OUT_FILENAME=[FILENAME,'_PhasevsR.dat'];
eval(string(strcat('save',{' '}, OUT_FILENAME,' C -ascii')));  
end

       % Coded in a quick and dirty way... For Info, Matteo x8498
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%