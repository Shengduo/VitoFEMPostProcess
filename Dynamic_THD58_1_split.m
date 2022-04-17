%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT RESULT FILES FROM VIC-2D
% 1. IMPORT MAT FILES
% 2. FILTER DISPLACEMENTS TWICE WITH NL-MEANS FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V.RUBINO

clc; clear all; close all

Test = 'THD58_1'; 

% LOADING
alpha = 29; % inclination angle of the interface
P = 14.3; % Applied load in MPa
sigma0 = P*(cosd(alpha))^2;
sigma0fp = P*(sind(alpha))^2;
tau0 = P*(sind(alpha))*(cosd(alpha));

% HIGH-SPEED CAMERA
camera_delay = 0;  % camera delay in us
nucl_delay = 8;  % nucleation delay in us (determined by trigger test)
rec_delay =  0; % us read from hpv-x after recording (< 0  when before)
frame_rate = 1;     % in Mfps 
inter = 1/frame_rate;        % interframe time in us
nof = 128; %number of files to be analyzed

time=inter*(0:nof-2); % if no delay time starts at 0us after trigger
time2=inter*(0:nof-1); 
time_delay = time + camera_delay + rec_delay - nucl_delay;  % to use with velocities
time_delay2 = time2 + camera_delay + rec_delay - nucl_delay; % to use with displacements 

% CORRELATION
Subset = 51;    
stepsize = 1;
Pix_10mm = 81;  %Number of pixels in a grid of 10 mm ;  144 pixels in 2 mm -> 72 pixels in 1 mm
pxsize=10/(Pix_10mm)*1000; % pixel size in um/pixel  (Pix_10mm pixels  in 10 mm)
%width = 18;     % Width of FOV
width = round(pxsize*400/1000);

pix_off = 6;       % pixels with no results
fov_dist = 39+17.5;   % distance in mm of fov from nucleation

pix_first = pix_off + 1;
% pix_last = size(FPdot_fil_twice,2) - pix_off;  
pix_last = 399 - pix_off;  
offset = pix_off * pxsize*10^-3;

% PHYSICAL PARAMETERS
cp = 2.61; % km/s or mm/us
cs = 1.29;
cR = 0.92*cs;
E = 5300; % MPa
vi = 0.35;

%% Files to read
[file_name dir_path]=uigetfile('*.mat','Pick the first file in the directory (.mat)');
F=dir([dir_path, '*.mat']);

% %%%% Check "Coordinates.m" script for the coordinates indices

for k=1:nof
    FILEtoREAD = [dir_path (F(k).name)];   % reads txt files and stores them in structure F
    results=load(FILEtoREAD);                % results is a structure
    x_up(:,:,k)=results.x;
    y_up(:,:,k)=results.y;
    x_low(:,:,k)=results.x_0;
    y_low(:,:,k)=results.y_0;
    sigma_up(:,:,k)=results.sigma;
    sigma_low(:,:,k)=results.sigma_0;
    EW_up(:,:,k)=pxsize*results.u;
    NS_up(:,:,k)=-pxsize*results.v;
    EW_low(:,:,k)=pxsize*results.u_0;
    NS_low(:,:,k)=-pxsize*results.v_0;
end

top_half = 1:size(EW_up,1);

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD/
mkdir ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)]);
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

%% Assembly unfiltered fields

kini=1;kfin=128;

Rows_abo=zeros(size(y_up(:,:,k)));
Rows_bel=zeros(size(y_low(:,:,k)));    

for k=kini:kfin
     EW_UP(:,:,k) = [EW_up(:,:,k);Rows_bel];
     NS_UP(:,:,k) = [NS_up(:,:,k);Rows_bel];
     EW_LOW(:,:,k) = [Rows_abo;EW_low(:,:,k)];
     NS_LOW(:,:,k) = [Rows_abo;NS_low(:,:,k)];
     EW(:,:,k) =  EW_UP(:,:,k) + EW_LOW(:,:,k);
     NS(:,:,k) =  NS_UP(:,:,k) + NS_LOW(:,:,k);
     FP(:,:,k) =  EW(:,:,k); % FP component in um
     FN(:,:,k) =  NS(:,:,k);% FN component in um 
end


for k=(kini+1):(kfin-1)
  FPdot(:,:,k) = (FP(:,:,k+1)-FP(:,:,k-1))/(2*inter);  % m/s=um/us
  FNdot(:,:,k) = (FN(:,:,k+1)-FN(:,:,k-1))/(2*inter);  % m/s=um/us
  Vel_mag(:,:,k) = sqrt(FPdot(:,:,k).^2 + FNdot(:,:,k).^2);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this script to have mm scale on the plots
xint = 10;      % interval of subdivisions in mm 
yint = xint;
xfirst = xint;  % first label to appear in mm
yfirst = yint;

xpix_int = xint/(pxsize*10^-3);   % interval of subdivisions in pixels 
ypix_int = yint/(pxsize*10^-3);
pix_first = xfirst/(pxsize*10^-3); % first pixel to be labelled
fov_width = (pxsize*10^-3)*400;    % fov width in mm 
fov_height = (pxsize*10^-3)*250;   % fov height in mm
numint = floor(fov_width/xint); 

xtick_loc = pix_first:xpix_int:400;  % location of labels in pixels 
xtick_lab = xfirst:xint:fov_width;   % label values

ytick_loc = pix_first:ypix_int:400;
ytick_lab = yfirst:yint:fov_height;

% Make movies of FP unfiltered

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
mkdir ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)]);
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])


writerObj = VideoWriter('FP_unfiltered.avi'); 
writerObj.FrameRate = 10;
open(writerObj);

% Set range
clims=[-30 30];

%figure
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 800]);
set(gcf,'color','w');

k1 = 70;%55;
k2 = 70;%80;
lim1 = 10;
lim2 = 30;
lim3 = 30;

for k=1:nof
    % Set range in intervals
    if k <= k1
        clims = [-lim1 lim1];
    elseif k > k1 & k <= k2 
        clims = [-lim2 lim2];
    elseif k > k2
        clims = [-lim3 lim3];
    end
    imagesc(FP(:,:,k),clims);
    axis equal
    axis tight
    set(gca,'fontsize',16)
    title(['Fault parallel displacement (\mum) - Frame = ' num2str(k)],'Fontsize',22)
    %title(['Fault parallel displacement (\mum) - time = ' num2str(time_delay2(k)) ' \mus'],'Fontsize',22)
    colorbar
    hcol=colorbar; 
    set(hcol,'fontsize',20)
    
    set(gca,'XTick',xtick_loc)
    set(gca,'XTickLabel',xtick_lab)
    set(gca,'fontsize',16)
    set(gca,'YTick',ytick_loc)
    set(gca,'YTickLabel',ytick_lab)
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Make movies of FN

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
mkdir ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)]);
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

writerObj = VideoWriter('FN_unfiltered.avi'); 
writerObj.FrameRate = 10;
open(writerObj);

% Set range
clims=[-5 5];

%figure
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 800]);
set(gcf,'color','w');

k1 = 40;
k2 = 130;
lim1 = 5;
lim2 = 10;
lim3 = 5;

for k=1:nof
%    Set range in intervals
    if k <= k1
        clims = [-lim1 lim1];
    elseif k > k1 & k <= k2 
        clims = [-lim2 lim2];
    elseif k > k2
        clims = [-lim3 lim3];
    end
    imagesc(FN(:,:,k),clims);
    axis equal
    axis tight
    set(gca,'fontsize',16)
    title(['Fault normal displacement (\mum) - Frame = ' num2str(k)],'Fontsize',22)
    colorbar
    hcol=colorbar; 
    set(hcol,'fontsize',20)
    
    set(gca,'XTick',xtick_loc)
    set(gca,'XTickLabel',xtick_lab)
    set(gca,'fontsize',16)
    set(gca,'YTick',ytick_loc)
    set(gca,'YTickLabel',ytick_lab)
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Make movies of FPdot
writerObj = VideoWriter('FPdot_unfiltered.avi'); 
writerObj.FrameRate = 10;
open(writerObj);

% Set range
clims=[-1 1];

k1 = 100;
k2 = 130;
lim1 = 1;
lim2 = 1;
lim3 = 1;

%figure
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 800]);
set(gcf,'color','w');

for k=1:nof-1  %90
    if k <= k1
        clims = [-lim1 lim1];
    elseif k > k1 & k <= k2 
        clims = [-lim2 lim2];
    elseif k > k2
        clims = [-lim3 lim3];
    end
    imagesc(FPdot(:,:,k),clims);
    axis equal
    axis tight
    set(gca,'fontsize',16)
    title(['Fault parallel velocity (m/s) - Frame = ' num2str(k)],'Fontsize',22)
    %title(['Fault parallel velocity (m/s) - time = ' num2str(time_delay(k))],'Fontsize',22)
    colorbar
    hcol=colorbar; 
    set(hcol,'fontsize',20)
        set(gca,'XTick',xtick_loc)
    set(gca,'XTickLabel',xtick_lab)
    set(gca,'fontsize',16)
    set(gca,'YTick',ytick_loc)
    set(gca,'YTickLabel',ytick_lab)
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Make movies of FNdot
writerObj = VideoWriter('FNdot_unfiltered.avi'); 
writerObj.FrameRate = 10;
open(writerObj);

% Set range
clims=[-1 1];

k1 = 100;
k2 = 130;
lim1 = 1;
lim2 = 1;
lim3 = 1;

%figure
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 800]);
set(gcf,'color','w');

for k=1:nof-1  %90
    if k <= k1
        clims = [-lim1 lim1];
    elseif k > k1 & k <= k2 
        clims = [-lim2 lim2];
    elseif k > k2
        clims = [-lim3 lim3];
    end
    imagesc(FNdot(:,:,k),clims);
    axis equal
    axis tight
    set(gca,'fontsize',16)
    title(['Fault parallel velocity (m/s) - Frame = ' num2str(k)],'Fontsize',22)
    %title(['Fault parallel velocity (m/s) - time = ' num2str(time_delay(k))],'Fontsize',22)
    colorbar
    hcol=colorbar; 
    set(hcol,'fontsize',20)
        set(gca,'XTick',xtick_loc)
    set(gca,'XTickLabel',xtick_lab)
    set(gca,'fontsize',16)
    set(gca,'YTick',ytick_loc)
    set(gca,'YTickLabel',ytick_lab)
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Specify filtering parameters
f=2; % patch size (3x3)
t=25; % search area dimension (21x21)
h=0.5;  % noise parameter

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/

kini = 1;  
kfin = 80;  %80

for k=kini:kfin   %40-75
  EW_fil_up(:,:,k)=NLmeans(EW_up(:,:,k),f,t,h);
  NS_fil_up(:,:,k)=NLmeans(NS_up(:,:,k),f,t,h);
  EW_fil_low(:,:,k)=NLmeans(EW_low(:,:,k),f,t,h);
  NS_fil_low(:,:,k)=NLmeans(NS_low(:,:,k),f,t,h);
  EW_fil2_up(:,:,k)=NLmeans(EW_fil_up(:,:,k),f,t,h);
  NS_fil2_up(:,:,k)=NLmeans(NS_fil_up(:,:,k),f,t,h);
  EW_fil2_low(:,:,k)=NLmeans(EW_fil_low(:,:,k),f,t,h);
  NS_fil2_low(:,:,k)=NLmeans(NS_fil_low(:,:,k),f,t,h);
end

figure('windowstyle','docked')
imagesc(EW_fil2_up(:,:,k),[-50 50]); axis equal; axis tight; colorbar

cd (['/Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD/' Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])
%save('THD7_all_split')
save('THD58_1_split_variables','EW_up','EW_low', 'NS_up', 'NS_low', ...
                            'EW_fil_up','NS_fil_up','EW_fil_low','NS_fil_low',...
                            'EW_fil2_up','NS_fil2_up','EW_fil2_low','NS_fil2_low')

% Continue filtering

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/

kini = 81;  %81
kfin = 128;  

for k=kini:kfin   %40-75
  EW_fil_up(:,:,k)=NLmeans(EW_up(:,:,k),f,t,h);
  NS_fil_up(:,:,k)=NLmeans(NS_up(:,:,k),f,t,h);
  EW_fil_low(:,:,k)=NLmeans(EW_low(:,:,k),f,t,h);
  NS_fil_low(:,:,k)=NLmeans(NS_low(:,:,k),f,t,h);
  EW_fil2_up(:,:,k)=NLmeans(EW_fil_up(:,:,k),f,t,h);
  NS_fil2_up(:,:,k)=NLmeans(NS_fil_up(:,:,k),f,t,h);
  EW_fil2_low(:,:,k)=NLmeans(EW_fil_low(:,:,k),f,t,h);
  NS_fil2_low(:,:,k)=NLmeans(NS_fil_low(:,:,k),f,t,h);
end

figure('windowstyle','docked')
imagesc(EW_fil2_up(:,:,k),[-50 50]); axis equal; axis tight; colorbar
title('EW fil2 up')

cd (['/Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD/' Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])
save('THD58_1_split_variables2','EW_up','EW_low', 'NS_up', 'NS_low', ...
                            'EW_fil_up','NS_fil_up','EW_fil_low','NS_fil_low',...
                            'EW_fil2_up','NS_fil2_up','EW_fil2_low','NS_fil2_low')
                        
%% Continue filtering 2

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/
f=2; % patch size (3x3)
t=25; % search area dimension (21x21)
h=0.5;  % noise parameter

kini = 106;  %54
kfin = 128;  %54

for k=kini:kfin   %40-75
  EW_fil_up(:,:,k)=NLmeans(EW_up(:,:,k),f,t,h);
  NS_fil_up(:,:,k)=NLmeans(NS_up(:,:,k),f,t,h);
  EW_fil_low(:,:,k)=NLmeans(EW_low(:,:,k),f,t,h);
  NS_fil_low(:,:,k)=NLmeans(NS_low(:,:,k),f,t,h);
  EW_fil2_up(:,:,k)=NLmeans(EW_fil_up(:,:,k),f,t,h);
  NS_fil2_up(:,:,k)=NLmeans(NS_fil_up(:,:,k),f,t,h);
  EW_fil2_low(:,:,k)=NLmeans(EW_fil_low(:,:,k),f,t,h);
  NS_fil2_low(:,:,k)=NLmeans(NS_fil_low(:,:,k),f,t,h);
end

figure('windowstyle','docked')
imagesc(EW_fil2_up(:,:,k),[-50 50]); axis equal; axis tight; colorbar
title(['EW fil2 up - ' num2str(k)])

cd (['/Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD/' Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])
%save({[num2str(Test) '_split_variables'],['EW_fil_up','NS_fil_up','EW_fil_low','NS_fil_low','EW_fil2_up','NS_fil2_up','EW_fil2_low','NS_fil2_low']})
save('THD58_1_split_variables2','EW_up','EW_low', 'NS_up', 'NS_low', ...
                            'EW_fil_up','NS_fil_up','EW_fil_low','NS_fil_low',...
                            'EW_fil2_up','NS_fil2_up','EW_fil2_low','NS_fil2_low')



%% Assembly filtered fields

kini=1;kfin=128;

Rows_abo=zeros(size(y_up(:,:,k)));
Rows_bel=zeros(size(y_low(:,:,k)));    

for k=kini:kfin
     EW_fil_UP(:,:,k) = [EW_fil_up(:,:,k);Rows_bel];
     NS_fil_UP(:,:,k) = [NS_fil_up(:,:,k);Rows_bel];
     EW_fil_LOW(:,:,k) = [Rows_abo;EW_fil_low(:,:,k)];
     NS_fil_LOW(:,:,k) = [Rows_abo;NS_fil_low(:,:,k)];    
          
     EW_fil2_UP(:,:,k) = [EW_fil2_up(:,:,k);Rows_bel];
     NS_fil2_UP(:,:,k) = [NS_fil2_up(:,:,k);Rows_bel];
     EW_fil2_LOW(:,:,k) = [Rows_abo;EW_fil2_low(:,:,k)];
     NS_fil2_LOW(:,:,k) = [Rows_abo;NS_fil2_low(:,:,k)];    
end

for k=kini:kfin
     Sigma_UP(:,:,k) = [sigma_up(:,:,k);Rows_bel];   
     Sigma_LOW(:,:,k) = [Rows_abo;sigma_low(:,:,k)];
end

Sigma = Sigma_UP + Sigma_LOW;

% Need also to add columns if EW and NS matrix have different number of columns

for k=kini:kfin   
    EW_fil(:,:,k) =  EW_fil_UP(:,:,k) + EW_fil_LOW(:,:,k);
    NS_fil(:,:,k) =  NS_fil_UP(:,:,k) + NS_fil_LOW(:,:,k);
    FP_fil(:,:,k) =  EW_fil(:,:,k); % FP component in um
    FN_fil(:,:,k) =  NS_fil(:,:,k);% FN component in um 
    
    EW_fil_twice(:,:,k) =  EW_fil2_UP(:,:,k) + EW_fil2_LOW(:,:,k);
    NS_fil_twice(:,:,k) =  NS_fil2_UP(:,:,k) + NS_fil2_LOW(:,:,k);
    FP_fil_twice(:,:,k) =  EW_fil_twice(:,:,k); % FP component in um
    FN_fil_twice(:,:,k) =  NS_fil_twice(:,:,k);% FN component in um 
end

% % Compute velocity field
% for k=(kini+1):(kfin-1)
%   FPdot_fil(:,:,k)=(FP_fil(:,:,k+1)-FP_fil(:,:,k-1))/(2*inter);  % m/s=um/us
%   FNdot_fil(:,:,k)=(FN_fil(:,:,k+1)-FN_fil(:,:,k-1))/(2*inter);  % m/s=um/us
%   Vel_mag_fil(:,:,k) = sqrt(FPdot_fil(:,:,k).^2 + FNdot_fil(:,:,k).^2);
% end


for k=(kini+1):(kfin-1)
  FPdot_fil_twice(:,:,k)=(FP_fil_twice(:,:,k+1)-FP_fil_twice(:,:,k-1))/(2*inter);  % m/s=um/us
  FNdot_fil_twice(:,:,k)=(FN_fil_twice(:,:,k+1)-FN_fil_twice(:,:,k-1))/(2*inter);  % m/s=um/us
  Vel_mag_fil_twice(:,:,k) = sqrt(FPdot_fil_twice(:,:,k).^2 + FNdot_fil_twice(:,:,k).^2);
end



%% Average top and bottom - SYMMETRY ADJUSTEMENT

% Check if top is larger than bottom half (rows_difference > 0)
% and shorten the larger domain in order to apply the symmetry adjustement
row_inter = size(EW_up,1);
top_rows_chk = 1:row_inter;
bot_rows_chk = (row_inter+1):size(FP_fil_twice,1);
rows_difference = size(FP_fil_twice(top_rows_chk,:,k),1) - size(FP_fil_twice(bot_rows_chk,:,k),1); % (top rows) - (bottom rows)
% rows_difference = size(top_rows_chk,2) - size(bot_rows_chk,2)


if rows_difference > 0
    top_rows = (rows_difference+1):row_inter;
    bot_rows = bot_rows_chk;
elseif rows_difference == 0
    top_rows = top_rows_chk;
    bot_rows = bot_rows_chk;
elseif rows_difference < 0
    top_rows = top_rows_chk;
    bot_rows = (row_inter+1):(size(FPdot_fil_twice,1)-(abs(rows_difference)));
end

FP_avg = zeros(size(FP_fil_twice));
FN_avg = zeros(size(FN_fil_twice));
FPdot_avg = zeros(size(FPdot_fil_twice));
FNdot_avg = zeros(size(FNdot_fil_twice));
Vel_mag_avg = zeros(size(FPdot_fil_twice));

clear counter
max_rows = size(top_rows,2);

for k=kini:kfin
   counter = 0;
   for i_top = top_rows(1):top_rows(end)
     % i_top is an index that runs from the image top to the interface
     % i_bot is an index that runs from the image bottom to the interface   
     index_bot = max_rows - counter;
     i_bot = bot_rows(index_bot);

     FP_avg(i_top,:,k) = (FP_fil_twice(i_top,:,k) - FP_fil_twice(i_bot,:,k))/2; %average of top with bottom changed of sign
     FP_avg(i_bot,:,k) = (FP_fil_twice(i_bot,:,k) - FP_fil_twice(i_top,:,k))/2; %average of top with bottom changed of sign

     FN_avg(i_top,:,k) = (FN_fil_twice(i_top,:,k) + FN_fil_twice(i_bot,:,k))/2; %average of top with bottom 
     FN_avg(i_bot,:,k) = (FN_fil_twice(i_bot,:,k) + FN_fil_twice(i_top,:,k))/2; %average of top with bottom 

     counter = counter +1;
   end
end

for k=(kini+1):(kfin-1)
  FPdot_avg(:,:,k)=(FP_avg(:,:,k+1)-FP_avg(:,:,k-1))/(2*inter);  % m/s=um/us
  FNdot_avg(:,:,k)=(FN_avg(:,:,k+1)-FN_avg(:,:,k-1))/(2*inter);  % m/s=um/us
  Vel_mag_avg(:,:,k) = sqrt(FPdot_avg(:,:,k).^2 + FNdot_avg(:,:,k).^2);
end


%% Clear and save data

clear ('EW_low','NS_up', 'NS_low','EW_fil_up','NS_fil_up','EW_fil_low','NS_fil_low','EW_fil2_up','NS_fil2_up','EW_fil2_low','NS_fil2_low');
clear EW_fil_UP NS_fil_UP EW_fil_LOW NS_fil_LOW EW_fil2_UP NS_fil2_UP EW_fil2_LOW NS_fil2_LOW Sigma_UP Sigma_LOW;
clear EW_fil NS_fil FP_fil FN_fil EW_fil_twice NS_fil_twice
clear FPdot_fil FNdot_fil Vel_mag_fil
clear Sigma_UP Sigma_LOW Sigma

save('THD58_1_Dynamic_sym.mat')




%% Check difference between filtered and unfiltered profiles

% Profile along the interface
% Select range along x
rangex = 14:393;

figure('windowstyle','docked')
plot(EW_up(size(EW_up,1),rangex,48),'k'); grid on; hold on
plot(EW_fil_up(size(EW_fil_up,1),rangex,48),'b')
plot(EW_fil2_up(size(EW_fil2_up,1),rangex,48),'r')
print -depsc EW_up_int_check.eps

%% Profile perpendicular to the interface
figure('windowstyle','docked')
plot(EW_up(20:size(EW_up,1),200,48),'k'); grid on; hold on
plot(EW_fil_up(20:size(EW_up,1),200,48),'b')
plot(EW_fil2_up(20:size(EW_up,1),200,48),'r')
print -depsc EW_up_perp_int_check.eps

%% GROUND MOTION ANALYSIS
% Movie of FPdot and FNdot vs. time on a path perpendicular to the fault

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

jset = 204;

writerObj = VideoWriter(['FPdot_along_i_jset_' num2str(jset) '.avi']); 
writerObj.FrameRate = 20;
open(writerObj);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 600]);
set(gcf,'color','w');

for i = (size(EW_up,1)+1):231   %size(EW_fil_twice,1)
    plot(time_delay,squeeze(FPdot_fil_twice(i,jset,:)),'b',time_delay,squeeze(FNdot_fil_twice(i,jset,:)),'r','Linewidth',2); 
    set(gca,'fontsize',17)
    xlabel('time (\mus)')
    ylabel('Particle velocity (m/s)')
    title(['jset = ' num2str(jset) ';  i = ' num2str(i)])
    grid on
    legend('FPdot','FNdot','Location','NorthWest')
    xlim([camera_delay max(time_delay)])
    ylim([-2 2])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Movie of FPdot and FNdot vs.time on a path along to the fault 

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

iset = size(EW_up,1);

l1 = -2;   
l2 = 2;

writerObj = VideoWriter(['FPdot_and_FNdot_along_fault.avi']); 
writerObj.FrameRate = 20;
open(writerObj);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 600]);
set(gcf,'color','w');

for j = (pix_first) : (pix_last)  
    plot(time_delay,squeeze(FPdot_fil_twice(iset,j,:)),'b',time_delay,squeeze(FNdot_fil_twice(iset,j,:)),'r','Linewidth',2); 
    set(gca,'fontsize',17)
    xlabel('time (\mus)')
    ylabel('Particle velocity (m/s)')
    title(['iset = ' num2str(iset) ';  j = ' num2str(j)])
    grid on
    legend('FPdot','FNdot','Location','NorthWest')
    xlim([min(time_delay) max(time_delay)])
    ylim([l1 l2])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Movie of FPdot and FNdot vs.time on a path along to the fault with arrival times; nucleaton delay included

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

iset = size(EW_up,1);

writerObj = VideoWriter(['FPdot_and_FNdot_along_fault_tarr.avi']); 
writerObj.FrameRate = 20;
open(writerObj);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 600]);
set(gcf,'color','w');

for j = (pix_first) : (pix_last)  
    hold off
    plot(time_delay,squeeze(FPdot_fil_twice(iset,j,:)),'b',time_delay,squeeze(FNdot_fil_twice(iset,j,:)),'r','Linewidth',2); 
        tp = ((j*pxsize/1000)+fov_dist)/cp;
        ts = ((j*pxsize/1000)+fov_dist)/cs;
        tR = ((j*pxsize/1000)+fov_dist)/cR;
        tsq2 = ((j*pxsize/1000)+fov_dist)/(sqrt(2)*cs);
        hold on
        plot([tp,tp],[l1 l2],'g')
        plot([ts,ts],[l1 l2],'g')
        plot([tR,tR],[l1 l2],'g')
        plot([tsq2,tsq2],[l1 l2],'g')
    set(gca,'fontsize',17)
    xlabel('time (\mus)')
    ylabel('Particle velocity (m/s)')
    title(['iset = ' num2str(iset) ';  j = ' num2str(j)])
    grid on
    legend('FPdot','FNdot','Location','NorthWest')
    xlim([min(time_delay) max(time_delay)])
    ylim([l1 l2])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Movie of FPdot and FNdot vs.time on a path along to the fault with arrival times no nucleaton delay included

cd /Users/vitorubino/Documents/MATLAB/VIC-2D/HPV-X/Figures_THD
cd ([Test '_Dynamic_' num2str(Subset) '_' num2str(stepsize)])

iset = size(EW_up,1);

writerObj = VideoWriter(['FPdot_and_FNdot_along_fault_tarr_nonucldelay.avi']); 
writerObj.FrameRate = 20;
open(writerObj);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 600]);
set(gcf,'color','w');

for j = (pix_first) : (pix_last)  
    hold off
    plot(time_delay+nucl_delay,squeeze(FPdot_fil_twice(iset,j,:)),'b',time_delay+nucl_delay,squeeze(FNdot_fil_twice(iset,j,:)),'r','Linewidth',2); 
        tp = ((j*pxsize/1000)+fov_dist)/cp;
        ts = ((j*pxsize/1000)+fov_dist)/cs;
        tR = ((j*pxsize/1000)+fov_dist)/cR;
        tsq2 = ((j*pxsize/1000)+fov_dist)/(sqrt(2)*cs);
        hold on
        plot([tp,tp],[l1 l2],'g')
        plot([ts,ts],[l1 l2],'g')
        plot([tR,tR],[l1 l2],'g')
        plot([tsq2,tsq2],[l1 l2],'g')
    set(gca,'fontsize',17)
    xlabel('time (\mus)')
    ylabel('Particle velocity (m/s)')
    title(['iset = ' num2str(iset) ';  j = ' num2str(j)])
    grid on
    legend('FPdot','FNdot','Location','NorthWest')
    xlim([min(time_delay+nucl_delay) max(time_delay+nucl_delay)])
    ylim([l1 l2])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

%% Movie of FPdot (above and below the fault) vs. time on a path perpendicular to the fault

jset = 204;

writerObj = VideoWriter(['FPdot_along_i_jset_' num2str(jset) '_sym.avi']); 
writerObj.FrameRate = 20;
open(writerObj);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1200, 600]);
set(gcf,'color','w');

hlim = size(EW_up,1);

for i = (hlim+1):231   %size(EW_fil_twice,1)
    plot(time_delay,squeeze(FPdot_fil_twice(i,jset,:)),'b',time_delay,squeeze(FPdot_fil_twice(2*hlim-i+1,jset,:)),'r','Linewidth',2); 
    set(gca,'fontsize',17)
    xlabel('time (\mus)')
    ylabel('Fault parallel velocity (m/s)')
    title(['jset = ' num2str(jset)])
    grid on
    legend(['i = ' num2str(i)])
    xlim([camera_delay max(time_delay)])
    ylim([-12 12])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj)

