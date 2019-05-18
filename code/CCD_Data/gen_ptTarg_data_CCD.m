%close all;
clear all;
addpath('helper_functions/');

%% INPUT PARAMETERS START HERE %%

disp_true = 0; % display true target locations in final images?
rand_targ = 0; % random target locations in imaging grid? if = 0 then linearly spaced
dimension = 2; % number of dimensions
spacing = 20; % spacing between points
data.N_targets = 24; % total number of points
data.iso = 1; % isotropic point scatterers?
data.num_ccd = 5; % number of points to change.
data.remove = 3; % number of points that are removed (the remaining change location)

% Define sensor parameters here
data.Fc = 9.6*1e9;            % Center freq (Hz)
data.xRes = .3;              % Cross Range Resolution at scene center (m)
data.rRes = .3;              % Range Resolution (m)
data.za = 5000;                % Altitude of sensor (m)
data.elev = (30)*pi/180;       % Depression angle at scene center (rad)
data.cent_angle = 0*180/pi;         % Center Angle (rad)
data.Np = 512;               % Number of pulses to image
data.K = 512;                % Number of frequency samples
data.Nx = 512;
data.Ny = 512;

% Define display parameters here
dyn_range = 50;              % dB of dynamic range to display
taylorFlag = 0;              % Add a Taylor window
slobe = 35;                  % Sidelobe level of Taylor window (dB)
nbar = 4;                    % n_bar for Taylor window
minDynRange = -30;

std_noise = .5; % standard deviation of AGWN (set to zero for no AGWN)


%% Gather imaging variables

% Calculate the constant range to sensor (m)
data.ra = data.za/sin(data.elev);

% Define speed of light (m/s)
c = 299792458;

% Calculate wavelength (m)
data.lambda = c / data.Fc;

% Calculate the frequency vector (Hz)
data.BW = c / (2*data.rRes);
data.freq = linspace(data.Fc-data.BW/2,data.Fc+data.BW/2,data.K);
data.spat_freq = 4*pi/c*data.freq'*cos(data.elev);

% Calculate length of synthetic aperture (rad)
data.intAngle = c / (2 * data.Fc * data.xRes);

% Initialize time vector (s)
t = linspace(-1,1,data.Np);

% Calculate circular flight path
data.azim = (t * data.intAngle / 2 + data.cent_angle)'; % (rad)
data.AntX = data.ra * cos(data.elev) * cos(data.azim);
data.AntY = data.ra * cos(data.elev) * sin(data.azim);
data.AntZ = data.za * ones(size(data.AntX));
data.R0 = sqrt(data.AntX.^2 + data.AntY.^2 + data.AntZ.^2);

% Calculate the frequency step size (Hz)
data.deltaF = diff(data.freq(1:2));

% Determine the average azimuth angle step size (radians)
data.deltaAz = abs(mean(diff(data.azim)));

% Determine the maximum scene size of the image (m)
data.maxWr = c/(2*data.deltaF) * cos(data.elev);
data.maxWx = c/(2*data.deltaAz*data.Fc);

data.pix_sizeR = data.maxWr/data.Ny;
data.pix_sizeX = data.maxWx/data.Nx;
data.x0 = 0;
data.y0 = 0;

%% Create scene

% reference scene
if dimension == 3
    [targPos,data] = gen_3D_pts(rand_targ,spacing,data);
elseif dimension == 2
    [targPos,data] = gen_2D_pts(rand_targ,spacing,data);
end

%% change detection scene
[targPos_CCD,data] = gen_CCD_pts(targPos,dimension,data);

%% Generate PHD

disp('Generating PHD...');

% Loop through each pulse to calculate the phase history
for ii = 1:length(data.azim)
    
    % Initialize the vector which contains the phase history for this pulse
    freqdata = zeros(1,data.K);
    % Loop through each target
    for kk = 1:data.N_targets
        % Calculate the differential range to the target (m)
        dR = sqrt((data.AntX(ii)-targPos(kk,1))^2+...
            (data.AntY(ii)-targPos(kk,2))^2+...
            (data.AntZ(ii)-targPos(kk,3))^2) - data.R0(ii);
        
        % Update the phase history for this pulse
        freqdata = freqdata + data.amp(kk) * exp(-1i*4*pi*dR/c*data.freq);
    end
    
    % Put the phase history into the data structure
    data.phdata(ii,:) = freqdata;
    
    
    freqdata_CCD = zeros(1,data.K);
    % Loop through each target
    for kk = 1:data.N_targets_CCD
        % Calculate the differential range to the target (m)
        dR = sqrt((data.AntX(ii)-targPos_CCD(kk,1))^2+...
            (data.AntY(ii)-targPos_CCD(kk,2))^2+...
            (data.AntZ(ii)-targPos_CCD(kk,3))^2) - data.R0(ii);
        
        % Update the phase history for this pulse
        freqdata_CCD = freqdata_CCD + data.amp_CCD(kk) * exp(-1i*4*pi*dR/c*data.freq);
    end
    
    % Put the phase history into the data structure
    data.phdata_CCD(ii,:) = freqdata_CCD;
    
    
    if mod(ii,100) == 0
        disp(['angle ', num2str(ii), ' out of ', num2str(data.Np)]);
    end
end

% windowing for side lobe reduction
if taylorFlag
    taperSlow = taylor(data.Np,slobe,nbar);
    taperFast = taylor(data.K,slobe,nbar);
    data.phdata = data.phdata .* (taperSlow*taperFast');
    data.phdata_CCD = data.phdata_CCD .* (taperSlow*taperFast');
end

% AGWN
noise = std_noise*(randn(data.K,data.Np) + 1i*randn(data.K,data.Np));
data.phdata = data.phdata' + noise;
data.phdata_reference = data.phdata; 

noise = std_noise*(randn(data.K,data.Np) + 1i*randn(data.K,data.Np));
data.phdata_CCD = data.phdata_CCD' + noise;


%% Form 2D image

data.AntAzim = 180/pi*data.azim';
data.AntElev = 180/pi*data.elev*ones(1,data.Np);
data.freq = col(data.freq);
data.Wx = data.maxWx;
data.Wy = data.maxWr;
data.Wz =0;
data.minF = min(data.freq)*ones(size(data.AntAzim));
data.Nfft = 2048;

data.x_vec = linspace(data.x0 - data.Wx/2, data.x0 + data.Wx/2, data.Nx);
data.y_vec = linspace(data.y0 - data.Wy/2, data.y0 + data.Wy/2, data.Ny);
[data.x_mat,data.y_mat] = meshgrid(data.x_vec,data.y_vec);
data.z_mat = zeros(size(data.x_mat));

% bp image reference 
fprintf('Forming reference image... \n');
data.phdata = data.phdata_reference; 
data_bp = bpBasicFarField(data);
img_bp_reference = data_bp.im_final;

% bp image CCD 
fprintf('Forming changed image... \n');
data.phdata = data.phdata_CCD; 
data_bp = bpBasicFarField(data);
img_bp_CCD = data_bp.im_final;

%% plot results

xtmp = linspace(-floor(data.maxWx/2),floor(data.maxWx/2),data.Nx);
ytmp = linspace(-floor(data.maxWr/2),floor(data.maxWr/2),data.Ny);

figure; 
subplot(1,2,1); 
imagesc(xtmp,ytmp,img_dB(img_bp_reference),[-50,0]);
hold on; colormap jet

if disp_true
    for ii = 1:data.N_targets
        plot(targPos(ii,1),targPos(ii,2),'wo');
    end
end
axis xy image
title('BP Image -- Reference')

subplot(1,2,2); 
imagesc(xtmp,ytmp,img_dB(img_bp_CCD),[-50,0]);
hold on; colormap jet

if disp_true
    for ii = 1:data.N_targets_CCD
        plot(targPos_CCD(ii,1),targPos_CCD(ii,2),'wo');
    end
end
axis xy image
title('BP Image -- CCD')

%% 3D BP to the actual target locations (just for fun)
if dimension == 3
    
    data.x_mat = targPos(:,1);
    data.y_mat = targPos(:,2);
    data.z_mat = targPos(:,3);
    
    data.phdata = data.phdata_reference;
    data_bp_3D = bpBasicFarField(data);
    
    img_ref = 20*log10(abs(data_bp_3D.im_final)./max(abs(data_bp_3D.im_final)));
    
    
    data.x_mat = targPos_CCD(:,1);
    data.y_mat = targPos_CCD(:,2);
    data.z_mat = targPos_CCD(:,3);
    
    data.phdata = data.phdata_CCD; 
    data_bp_3D = bpBasicFarField(data);
    
    img_CCD = 20*log10(abs(data_bp_3D.im_final)./max(abs(data_bp_3D.im_final)));
    
    
    figure;
    subplot(1,2,1); 
    scatter3(targPos(:,1),targPos(:,2),targPos(:,3),50,...
        img_ref,'*');
    axis equal
    axis tight
    title('Reference')
    
    subplot(1,2,2); 
    scatter3(targPos_CCD(:,1),targPos_CCD(:,2),targPos_CCD(:,3),50,...
        img_CCD,'*');
    axis equal
    axis tight
    title('Changed')
    
end


