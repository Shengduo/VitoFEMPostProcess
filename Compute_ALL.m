%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE SLIP AND SLIP RATE
% COMPUTE STRAINS, STRESSES AND FRICTION
% PRODUCE FRICTION MOVIES

warning off;

% Compute slip and slip rate from filtered fields
row_inter = size(EW_up,1);  %121;
jstart = 10;
jend = 387; %390;
col_range = jstart:jend; % select x axis of sxy vs. x plot

Slip_avg = squeeze(FP_avg(row_inter+1,col_range,:) - FP_avg(row_inter,col_range,:));
Slip_rate_avg = squeeze(FPdot_avg(row_inter+1,col_range,:) - FPdot_avg(row_inter,col_range,:));


% Calculate strains after filtering twice - AVG
[M,N]=size(FP_avg(:,:,1));

hs = 1; % step to compute strains
bw = size(EW_up,1);   % vertical size of top half
bl = size(FP_avg,1);  % vertical size of full matrix

%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%
% Removing this comment will result in overwritng 
% the symetry adjusted variables
% Use only to deliberately compute raw quantities
% that are not symmetry adjusted
% FP_avg = FP_fil_twice;
% FN_avg = FN_fil_twice;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eFP_avg = zeros(size(FP_avg,1)-2*hs,size(FP_avg,2)-2*hs,size(FP_avg,3));
eFN_avg = zeros(size(FP_avg,1)-2*hs,size(FP_avg,2)-2*hs,size(FP_avg,3));
ePN_avg = zeros(size(FP_avg,1)-2*hs,size(FP_avg,2)-2*hs,size(FP_avg,3));


for k=kini:kfin
  for j = (1+2*hs) : (N-2*hs)
    for i = (1+2*hs) : (M-2*hs)  %bw
       if i >= 1 & i <= (bw-2*hs) 
           %Use central difference approximation                      
           eFP_avg(i,j,k)=(FP_avg(i,j+hs,k)-FP_avg(i,j-hs,k))/(2*hs*stepsize)/pxsize;   % x positive rightward
           eFN_avg(i,j,k)=-(FN_avg(i+hs,j,k)-FN_avg(i-hs,j,k))/(2*hs*stepsize)/pxsize;  % y positive upward
           ePN_avg(i,j,k)=0.5*(-(FP_avg(i+hs,j,k)-FP_avg(i-hs,j,k))+(FN_avg(i,j+hs,k)-FN_avg(i,j-hs,k)))/(2*hs*stepsize*pxsize);
           %ePN_avg_c1(i,j,k) = -(FP_avg(i+hs,j,k)-FP_avg(i-hs,j,k))/(2*hs*stepsize*pxsize); 
           %ePN_avg_c2(i,j,k) = (FN_avg(i,j+hs,k)-FN_avg(i,j-hs,k))/(2*hs*stepsize*pxsize); 
           
       elseif i > (bw-2*hs) & i <= bw   
           % Use second backward difference
           eFP_avg(i,j,k)=(FP_avg(i,j-2*hs,k)-4*FP_avg(i,j-hs,k)+3*FP_avg(i,j,k))/(2*hs*stepsize)/pxsize;   % x positive rightward
           eFN_avg(i,j,k)=-(FN_avg(i-2*hs,j,k)-4*FN_avg(i-hs,j,k)+3*FN_avg(i,j,k))/(2*hs*stepsize)/pxsize;  % y positive upward
           ePN_avg(i,j,k)=0.5*(-(FP_avg(i-2*hs,j,k)-4*FP_avg(i-hs,j,k)+3*FP_avg(i,j,k))+(FN_avg(i,j-2*hs,k)-4*FN_avg(i,j-hs,k)+3*FN_avg(i,j,k)))/(2*hs*stepsize*pxsize);
           %ePN_avg_c1(i,j,k) = -(FP_avg(i-2*hs,j,k)-4*FP_avg(i-hs,j,k)+3*FP_avg(i,j,k))/(2*hs*stepsize*pxsize);
           %ePN_avg_c2(i,j,k) = (FN_avg(i,j-2*hs,k)-4*FN_avg(i,j-hs,k)+3*FN_avg(i,j,k))/(2*hs*stepsize*pxsize);
           
       elseif i > bw & i <= (bw+2*hs)
           % Use second forward difference
           eFP_avg(i,j,k)=(-FP_avg(i,j+2*hs,k)+4*FP_avg(i,j+hs,k)-3*FP_avg(i,j,k))/(2*hs*stepsize)/pxsize;   % x positive rightward
           eFN_avg(i,j,k)=-(-FN_avg(i+2*hs,j,k)+4*FN_avg(i+hs,j,k)-3*FN_avg(i,j,k))/(2*hs*stepsize)/pxsize;  % y positive upward
           ePN_avg(i,j,k)=0.5*(-(-FP_avg(i+2*hs,j,k)+4*FP_avg(i+hs,j,k)-3*FP_avg(i,j,k))+(-FN_avg(i,j+2*hs,k)+4*FN_avg(i,j+hs,k)-3*FN_avg(i,j,k)))/(2*hs*stepsize*pxsize);
           %ePN_avg_c1(i,j,k) = -(-FP_avg(i+2*hs,j,k)+4*FP_avg(i+hs,j,k)-3*FP_avg(i,j,k))/(2*hs*stepsize*pxsize);
           %ePN_avg_c2(i,j,k) = (-FN_avg(i,j+2*hs,k)+4*FN_avg(i,j+hs,k)-3*FN_avg(i,j,k))/(2*hs*stepsize*pxsize);
           
       elseif i > (bw+2*hs) & i <= bl
           % Use central difference
           eFP_avg(i,j,k)=(FP_avg(i,j+hs,k)-FP_avg(i,j-hs,k))/(2*hs*stepsize)/pxsize;   % x positive rightward
           eFN_avg(i,j,k)=-(FN_avg(i+hs,j,k)-FN_avg(i-hs,j,k))/(2*hs*stepsize)/pxsize;  % y positive upward
           ePN_avg(i,j,k)=0.5*(-(FP_avg(i+hs,j,k)-FP_avg(i-hs,j,k))+(FN_avg(i,j+hs,k)-FN_avg(i,j-hs,k)))/(2*hs*stepsize*pxsize);
           %ePN_avg_c1(i,j,k) = -(FP_avg(i+hs,j,k)-FP_avg(i-hs,j,k))/(2*hs*stepsize*pxsize); 
           %ePN_avg_c2(i,j,k) = (FN_avg(i,j+hs,k)-FN_avg(i,j-hs,k))/(2*hs*stepsize*pxsize); 
          
       end
    end
  end
end




% COMPUTATION OF STRAIN RATES AVG

e_33_avg = vi/(1-vi)*(eFP_avg + eFN_avg);  % Assuming plane stress

for k=(kini+1):(kfin-1)
  eFPdot_avg(:,:,k) = (eFP_avg(:,:,k+1)-eFP_avg(:,:,k-1))/(2*inter*10^-6); %*10^-3;  % s^-1; 10^-6 to convert inter from us to s
  eFNdot_avg(:,:,k) = (eFN_avg(:,:,k+1)-eFN_avg(:,:,k-1))/(2*inter*10^-6); %*10^-3;  % 10^-3 is to plot strain rate in units (x 10^3) 
  ePNdot_avg(:,:,k) = (ePN_avg(:,:,k+1)-ePN_avg(:,:,k-1))/(2*inter*10^-6); %*10^-3;
  e_33dot_avg(:,:,k) = (e_33_avg(:,:,k+1)-e_33_avg(:,:,k-1))/(2*inter*10^-6);
end

Volum_strain_rate = eFPdot_avg + eFNdot_avg + e_33dot_avg;
edot_mag = sqrt(eFPdot_avg.^2 + eFNdot_avg.^2 + e_33dot_avg.^2 + ePNdot_avg.^2);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stresses - AVG
Ed = 5300;

sFP_avg = Ed/(1-vi^2)*(eFP_avg+vi*eFN_avg);
sFN_avg = Ed/(1-vi^2)*(vi*eFP_avg+eFN_avg);
sPN_avg = Ed/(1-vi^2)*(1-vi)*ePN_avg;    
tau_max_avg = sqrt(((sFP_avg - sFN_avg)/2).^2+sPN_avg.^2);   

% Compute stress changes along the interface
% Note that actual stresses and actual stresses along interface are computed below

row_inter_gap = size(EW_up,1);
%col_range = 10:390;

% Compute actual stress above interface
sPN_line_above_avg = squeeze(sPN_avg(row_inter_gap,col_range,:));
sFN_line_above_avg = squeeze(sFN_avg(row_inter_gap,col_range,:));
sFP_line_above_avg = squeeze(sFP_avg(row_inter_gap,col_range,:));

% Compute actual stress below interface
sPN_line_below_avg = squeeze(sPN_avg(row_inter_gap+1,col_range,:));
sFN_line_below_avg = squeeze(sFN_avg(row_inter_gap+1,col_range,:));
sFP_line_below_avg = squeeze(sFP_avg(row_inter_gap+1,col_range,:));

% Compute actual stress as average of above and below interface
sPN_line_AVG = (sPN_line_above_avg + sPN_line_below_avg)/2;
sFN_line_AVG = (sFN_line_above_avg + sFN_line_below_avg)/2;
sFP_line_AVG = (sFP_line_above_avg + sFP_line_below_avg)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute out-of-plane displacement, strains and stresses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_thick = 9.8; % thickness in mm 

eps_33_avg = -vi/(2*Ed)*(sFP_avg + sFN_avg);
u3_avg = -vi*h_thick/(2*Ed)*(sFP_avg + sFN_avg)*10^3;  % displacement expressed in microns

u3_avg_1 = -vi*h_thick/(2*Ed)*(sFP_avg)*10^3; 
u3_avg_2 = -vi*h_thick/(2*Ed)*(sFN_avg)*10^3; 

% u3dot_avg = zeros(size(FPdot_avg));
% u3dot_avg_1 = zeros(size(FPdot_avg));
% u3dot_avg_2 = zeros(size(FPdot_avg));

% Compute out-of-plane velocity
for k=(kini+1):(kfin-1)
   u3dot_avg(:,:,k) = (u3_avg(:,:,k+1) - u3_avg(:,:,k-1))/(2*inter);  % m/s=um/us
   u3dot_avg_1(:,:,k) = (u3_avg_1(:,:,k+1) - u3_avg_1(:,:,k-1))/(2*inter);
   u3dot_avg_2(:,:,k) = (u3_avg_2(:,:,k+1) - u3_avg_2(:,:,k-1))/(2*inter);
end


% Compute velocity magnitude including out-of-plane component
Vel_mag_avg_all = zeros(size(u3dot_avg));
for k=(kini+1):(kfin-1)
  Vel_mag_avg_all(:,:,k) = sqrt(FPdot_avg(2:end-1,2:end-1,k).^2 + FNdot_avg(2:end-1,2:end-1,k).^2 +  u3dot_avg(:,:,k).^2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual stresses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P
sigma0
sigma0fp
tau0

% Avg
sPN_actual_avg = tau0*ones(size(sPN_avg)) - sPN_avg;
sFN_actual_avg = sigma0*ones(size(sFN_avg)) - sFN_avg;
sFP_actual_avg = sigma0fp*ones(size(sFP_avg)) - sFP_avg;


% Compute actual stresses along the interface

row_inter_gap = size(EW_up,1);
% col_range = 10:390;

% Compute actual stress above interface
sPN_actual_line_above_avg = squeeze(sPN_actual_avg(row_inter_gap,col_range,:));
sFN_actual_line_above_avg = squeeze(sFN_actual_avg(row_inter_gap,col_range,:));
sFP_actual_line_above_avg = squeeze(sFP_actual_avg(row_inter_gap,col_range,:));

% Compute actual stress below interface
sPN_actual_line_below_avg = squeeze(sPN_actual_avg(row_inter_gap+1,col_range,:));
sFN_actual_line_below_avg = squeeze(sFN_actual_avg(row_inter_gap+1,col_range,:));
sFP_actual_line_below_avg = squeeze(sFP_actual_avg(row_inter_gap+1,col_range,:));

% Compute actual stress as average of above and below interface
sPN_actual_line_AVG = (sPN_actual_line_above_avg + sPN_actual_line_below_avg)/2;
sFN_actual_line_AVG = (sFN_actual_line_above_avg + sFN_actual_line_below_avg)/2;
sFP_actual_line_AVG = (sFP_actual_line_above_avg + sFP_actual_line_below_avg)/2;

% Compute friction coefficient
friction_ABO = sPN_actual_line_above_avg./sFN_actual_line_above_avg;
friction_BEL = sPN_actual_line_below_avg./sFN_actual_line_below_avg;
friction_AVG = sPN_actual_line_AVG./sFN_actual_line_AVG;

%% Compute stresses and friction away from the interface

half_subset = (Subset-1)/2;
i_above_awy = (row_inter-half_subset):row_inter;
i_below_awy = (row_inter+1):(row_inter+half_subset+1);

% Compute stresses away from the interface - above
sPN_actual_line_above_avg_awy = squeeze(sPN_actual_avg(i_above_awy,col_range,:));
sFN_actual_line_above_avg_awy = squeeze(sFN_actual_avg(i_above_awy,col_range,:));
sFP_actual_line_above_avg_awy = squeeze(sFP_actual_avg(i_above_awy,col_range,:));

% Compute stresses away from the interface - below
sPN_actual_line_below_avg_awy = squeeze(sPN_actual_avg(i_below_awy,col_range,:));
sFN_actual_line_below_avg_awy = squeeze(sFN_actual_avg(i_below_awy,col_range,:));
sFP_actual_line_below_avg_awy = squeeze(sFP_actual_avg(i_below_awy,col_range,:));

sPN_actual_line_AVG_awy = zeros(size(sPN_actual_line_above_avg_awy));
sFN_actual_line_AVG_awy = zeros(size(sPN_actual_line_above_avg_awy));
sFP_actual_line_AVG_awy = zeros(size(sPN_actual_line_above_avg_awy));

% Compute actual stress as average of above and below interface
for ii = 1:size(i_above_awy,2)
    sPN_actual_line_AVG_awy(ii,:,:) = (sPN_actual_line_above_avg_awy(ii,:,:) + sPN_actual_line_below_avg_awy(size(i_above_awy,2)+1-ii,:,:))/2;
    sFN_actual_line_AVG_awy(ii,:,:) = (sFN_actual_line_above_avg_awy(ii,:,:) + sFN_actual_line_below_avg_awy(size(i_above_awy,2)+1-ii,:,:))/2;
    sFP_actual_line_AVG_awy(ii,:,:) = (sFP_actual_line_above_avg_awy(ii,:,:) + sFP_actual_line_below_avg_awy(size(i_above_awy,2)+1-ii,:,:))/2;
end


% Compute friction coefficient
friction_ABO_awy = sPN_actual_line_above_avg_awy./sFN_actual_line_above_avg_awy;
friction_BEL_awy = sPN_actual_line_below_avg_awy./sFN_actual_line_below_avg_awy;
friction_AVG_awy = sPN_actual_line_AVG_awy./sFN_actual_line_AVG_awy;


%% Compute friction at a location away from the interface


% % Set far-field position
% % Use for small FOV
% % x2_farfield = 4; % mm
% % % Use for intermediate FOV
% x2_farfield = 10; % mm
% % % Use for large FOV
% % x2_farfield = 40; % mm
% px_farfield = x2_farfield/pxsize*10^3;
% i_far_above = row_inter - px_farfield;
% i_far_below = (row_inter+1) + px_farfield;
% 
% 
% % Compute actual stress above interface
% sPN_actual_line_above_avg_FF = squeeze(sPN_actual_avg(i_far_above,col_range,:));
% sFN_actual_line_above_avg_FF = squeeze(sFN_actual_avg(i_far_above,col_range,:));
% sFP_actual_line_above_avg_FF = squeeze(sFP_actual_avg(i_far_above,col_range,:));
% 
% % Compute actual stress below interface
% sPN_actual_line_below_avg_FF = squeeze(sPN_actual_avg(i_far_below,col_range,:));
% sFN_actual_line_below_avg_FF = squeeze(sFN_actual_avg(i_far_below,col_range,:));
% sFP_actual_line_below_avg_FF = squeeze(sFP_actual_avg(i_far_below,col_range,:));
% 
% % Compute actual stress as average of above and below interface
% sPN_actual_line_AVG_FF = (sPN_actual_line_above_avg_FF + sPN_actual_line_below_avg_FF)/2;
% sFN_actual_line_AVG_FF = (sFN_actual_line_above_avg_FF + sFN_actual_line_below_avg_FF)/2;
% sFP_actual_line_AVG_FF = (sFP_actual_line_above_avg_FF + sFP_actual_line_below_avg_FF)/2;
% 
% 
% % Compute friction coefficient
% friction_AVG_FF = sPN_actual_line_AVG_FF./sFN_actual_line_AVG_FF;
% Slip_avg_FF = squeeze(FP_avg(i_far_below,col_range,:) - FP_avg(i_far_above,col_range,:));
% Slip_rate_avg_FF = squeeze(FPdot_avg(i_far_below,col_range,:) - FPdot_avg(i_far_above,col_range,:));


%% Friction averaged along the interface

% friction_AVG_A = size(Slip_rate_avg_FF,2);
% Slip_avg_A = size(Slip_rate_avg_FF,2);
% Slip_rate_avg_A = size(Slip_rate_avg_FF,2);
% 
% friction_AVG_FFA = size(Slip_rate_avg_FF,2);
% Slip_avg_FFA = size(Slip_rate_avg_FF,2);
% Slip_rate_avg_FFA = size(Slip_rate_avg_FF,2);
% 
% for k = 1 : size(Slip_rate_avg_FF,2)
%    
%    % Average along the fault
%    friction_AVG_A(k) = mean(friction_AVG(:,k));
%    Slip_avg_A(k) = mean(Slip_avg(:,k));
%    Slip_rate_avg_A(k) = mean(Slip_rate_avg(:,k));
%    
%    % Average in the far-field
%    friction_AVG_FFA(k) = mean(friction_AVG_FF(:,k));
%    Slip_avg_FFA(k) = mean(Slip_avg_FF(:,k));
%    Slip_rate_avg_FFA(k) = mean(Slip_rate_avg_FF(:,k));
%    
% end
% 
% % Average in the far field
% friction_AVG_FFA_mean = mean(friction_AVG_FFA);
% Slip_avg_FFA_mean = mean(Slip_avg_FFA);
% Slip_rate_avg_FFA_mean = mean(Slip_rate_avg_FFA);


