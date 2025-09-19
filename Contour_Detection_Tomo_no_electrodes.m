%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CONTOUR DETECTION TOMOGRAPHY                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\'); % Working directory = directory where images are stored

TruncDir = [WorkDir '\1_Truncation']; % Directory for truncated images
CorrecDir = [WorkDir '\2.1_Correc']; % Directory for laser corrected images
AdjustDir = [WorkDir '\2.2_Adjust']; % Directory for contrast adjusted images
BinDir = [WorkDir '\2.3_Binary']; % Directory for binary images
mkdir(BinDir)
ContourDir = [WorkDir '\3_Contour']; % Directory for contour data
mkdir(ContourDir)


%% PARAMETERS

load([WorkDir '\Parameters']) % Load processing parameters
load([LaserDir '\Laser_Profile.mat']) % Load the laser intensity profile

HIT_rad = 20; % HIT zone radius [mm]
HIT_rad_pxl = floor(HIT_rad/Magn); % HIT zone radius [pxl]
HIT = true; % True by default, control if a flame is within the HIT zone

nbIm_start = 1; % First image to be processed for contour detection

% Erosion parameters
smfilt = 7; % Size of the median filter
erosion_size = 2;
se = strel('disk',erosion_size,6); % Erosion structure

ScaleFilter = 0.4; % Scale of filtering the pixelisation noise (mm)

%% BACKGROUND SUBSTRACTION

List_RawIm = dir([WorkDir '\*.tif']); % Gather all raw images in the directory
BACKGR = exist([WorkDir '\background.tif'], 'file'); % Checks if a backgorund image already exists
if BACKGR==0
    disp('#        NO BACKGROUND IMAGE FOUND        #')
    return
else
    List_RawIm = List_RawIm(2:end); % Removing background.tif from the list 
    Background = imread([WorkDir '\background.tif']); % Reading the background image
end

nbIm = size(List_RawIm,1); % Number of raw images
Name_RawIm = sortrows(char(List_RawIm.name)); % Name of the different images

TRUNC = exist(TruncDir,"dir"); % Check if truncated image directory exist
% if not, create it and truncate all images
if TRUNC == 0
    disp('### TRUNCATION ###')
    mkdir(TruncDir)
    tic 
    for i = 1:nbIm
        disp(['Truncation ' Name_RawIm(i,:) '...']);
        RawIm = imread([WorkDir '\' Name_RawIm(i,:)]);
        TruncIm = (imsubtract(RawIm,Background)); % Substracting the background 
        imwrite(TruncIm,[TruncDir '\Trunc_' Name_RawIm(i,end-8:end-4) '.tif']); % Saving the truncated image
    end
    disp('### TRUNCATION DONE ###')
    toc
end

List_TruncIm = dir([TruncDir '\*.tif']); % Gather all truncated images
Name_TruncIm = sortrows(char(List_TruncIm.name)); % Name of all gathered images


%% LASER SHEET CORRECTION

CORREC = exist(CorrecDir,'dir'); % Check if laser corrected images directory exist
% if not, create it and correct all truncated images
if CORREC == 0
    disp('### LASER CORRECTION ###')
    mkdir(CorrecDir)
    tic
    for i = 1:nbIm
        disp(['Correction ' Name_TruncIm(i,:) '...'])
        CorrecIm = im2double(imread([TruncDir '\' Name_TruncIm(i,:)]));
        for coll = 1:sizex
            CorrecIm(:,coll) = CorrecIm(:,coll)./(laserProfile_filt/max(laserProfile_filt)); % Normalisation to take account of the non-homogeneity of the laser sheet
        end
        imwrite(CorrecIm,[CorrecDir '\Correc_' Name_TruncIm(i,end-8:end-4) '.tif']); % Saving the corrected image
    end
    disp('### LASER CORRECTION DONE ###')
    toc
end

List_CorrecIm = dir([CorrecDir '\*.tif']); % Gather all laser corrected images
Name_CorrecIm = sortrows(char(List_CorrecIm.name)); % Name of all the laser corrected images


%% CONTRAST ADJUSTMENT

ADJUST = exist(AdjustDir, 'dir'); % Check if constrast adjusted images directory exist
% if not, create it and contrast adjust all laser corrected images
if ADJUST == 0
    disp('### CONTRAST ADJUSTMENT ###')
    mkdir(AdjustDir)
    tic
    for i = 1:nbIm
        disp(['Adjustment ' Name_CorrecIm(i,:) '...'])
        AdjustIm = im2double(imread([CorrecDir '\' Name_CorrecIm(i,:)]));
        AdjustIm = imadjust(AdjustIm); % Contrast adjustment of the image
        imwrite(AdjustIm,[AdjustDir '\Adjust_' Name_CorrecIm(i,end-8:end-4) '.tif']); % Saving the adjusted image
    end
    disp('### CONTRAST ADJUSTMENT DONE ###')
    toc
end

List_AdjustIm = dir([AdjustDir '\*.tif']); % Gather all the contrast adjusted images
Name_AdjustIm = sortrows(char(List_AdjustIm.name)); % Name of all the contrast adjusted images

%% MASKS

load viewport_mask.mat 
figure,imshow(imread([WorkDir '\' Name_RawIm(end,:)])) % Show the last raw image 
h1 = drawpolygon(gca, 'Position', viewport_pos); % Draw the mask for the viewport
wait(h1); 
disp('#    LASER SHEET MASK OK     #')

Mask_viewport = createMask(h1); % Create a binary mask for the viewport

viewport_pos = h1.Position; % Extracting position
save('viewport_mask.mat','viewport_pos')

close all

fig = figure('Visible','off');
imshow(imread([WorkDir '\' Name_RawIm(end,:)]))
h2 = drawcircle(gca,'Center',HIT_center,'Radius', HIT_rad_pxl); % Draw a circle corresponding to the HIT zone
Mask_HIT = createMask(h2); % Create a binary mask for the HIT zone
cont_HIT = imcontour(Mask_HIT,1); % Contour of the HIT zone

close all

%% BINARY THRESHOLD 

threshold_test = true; % Controls the loop
coeff = 1.0; % Coefficient to adjust the threshold

RawIm_first = imread([WorkDir '\' Name_RawIm(nbIm_start,:)]);
RawIm_last = imread([WorkDir '\' Name_RawIm(end,:)]);
AdjustIm_first = im2double(imread([AdjustDir '\' Name_AdjustIm(nbIm_start,:)]));
AdjustIm_last = im2double(imread([AdjustDir '\' Name_AdjustIm(end,:)]));

thresholdtemp = graythresh(AdjustIm_first); % temporary threshold used for the loop
fprintf('Binary threshold = %0.03f', thresholdtemp); % Affiche la valeur de seuil trouvée

while threshold_test
    thresholdtemp = thresholdtemp*coeff; % Updating the threshold

    % Full binarisation process for the first and last images
    BinIm_first = imbinarize(AdjustIm_first, thresholdtemp);
    BinIm_first = iminv(BinIm_first); 
    BinIm_first = BinIm_first.*Mask_viewport;
    % BinIm_first = imerode(BinIm_first,se); 
    % BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = Fct_Struct_Max(BinIm_first); 
    BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]);
    % BinIm_first = imdilate(BinIm_first,se);
    % BinIm_first = imdilate(BinIm_first,se);
    % BinIm_first = imdilate(BinIm_first,se);
    % BinIm_first = imdilate(BinIm_first,se);
    % BinIm_first = imdilate(BinIm_first,se);
    % BinIm_first = imerode(BinIm_first,se);
    % BinIm_first = imerode(BinIm_first,se);
    % BinIm_first = imerode(BinIm_first,se);
    % BinIm_first = imerode(BinIm_first,se);
    % BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imfill(BinIm_first,'holes'); 
    BinIm_first = Fct_Struct_Max(BinIm_first); 
    BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]); 

    BinIm_last = imbinarize(AdjustIm_last, thresholdtemp);
    BinIm_last = iminv(BinIm_last); 
    BinIm_last = BinIm_last.*Mask_viewport;
    % BinIm_last = imerode(BinIm_last,se); 
    % BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = Fct_Struct_Max(BinIm_last); 
    BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]);
    % BinIm_last = imdilate(BinIm_last,se);
    % BinIm_last = imdilate(BinIm_last,se);
    % BinIm_last = imdilate(BinIm_last,se);
    % BinIm_last = imdilate(BinIm_last,se);
    % BinIm_last = imdilate(BinIm_last,se);
    % BinIm_last = imerode(BinIm_last,se);
    % BinIm_last = imerode(BinIm_last,se);
    % BinIm_last = imerode(BinIm_last,se);
    % BinIm_last = imerode(BinIm_last,se);
    % BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imfill(BinIm_last,'holes'); 
    BinIm_last = Fct_Struct_Max(BinIm_last); 
    BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]); 

    % Contour detection for the first and last image to check if binarisation threshold is good
    fig_cont_first = figure('Visible','off');
    cont_first = imcontour(BinIm_first,1); % Contour detection
    close(fig_cont_first)

    fig_cont_last = figure('Visible','off');
    cont_last = imcontour(BinIm_last,1); % Contour detection
    close(fig_cont_last)

    x_cont_first = cont_first(1,2:end); % Contour coordinates
    y_cont_first = cont_first(2,2:end);

    x_cont_last = cont_last(1,2:end); % Contour coordinates
    y_cont_last = cont_last(2,2:end);

    % Check if contour is closed
    if x_cont_first(end)==x_cont_first(1) && y_cont_first(end) == y_cont_first(1)
        if x_cont_first(1)== x_cont_first(end) % If contour extremities are aligned
            x_cont_first = [x_cont_first x_cont_first(1)]; % closing manually the contour
            y_cont_first = [y_cont_first y_cont_first(1)];
        else
            slope = (y_cont_first(end)-y_cont_first(1))/(x_cont_first(end)-x_cont_first(1)); % Closing the contour using linear segment
            y_intersect = y_cont_first(1)-slope*x_cont_first(1);
            x_line = x_cont_first(1)+0.5:0.5:x_cont_first(end)-0.5;
            y_line = slope.*x_line + y_intersect;
            x_cont_first = [x_cont_first fliplr(x_line) x_cont_first(1)];
            y_cont_first = [y_cont_first fliplr(y_line) y_cont_first(1)];
        end
    end

    if x_cont_last(end)==x_cont_last(1) && y_cont_last(end) == y_cont_last(1)
        if x_cont_last(1)== x_cont_last(end) % If contour extremities are aligned
            x_cont_last = [x_cont_last x_cont_last(1)]; % closing manually the contour
            y_cont_last = [y_cont_last y_cont_last(1)];
        else
            slope = (y_cont_last(end)-y_cont_last(1))/(x_cont_last(end)-x_cont_last(1)); % Closing the contour using linear segment
            y_intersect = y_cont_last(1)-slope*x_cont_last(1);
            x_line = x_cont_last(1)+0.5:0.5:x_cont_last(end)-0.5;
            y_line = slope.*x_line + y_intersect;
            x_cont_last = [x_cont_last fliplr(x_line) x_cont_last(1)];
            y_cont_last = [y_cont_last fliplr(y_line) y_cont_last(1)];
        end
    end

    contour_first = x_cont_first -1i * y_cont_first; % Reconstructing the contour geometry
    contour_last = x_cont_last -1i * y_cont_last;

    % Contour filtering
    contour_filt_first = Fct_Contour_Filter(contour_first, ScaleFilter, Magn);
    contour_filt_last = Fct_Contour_Filter(contour_last, ScaleFilter, Magn);

    % Contour interpolation
    % Compute the distance between two successive points 
    % If interpolation is good the distance should be around 0.5 
    % Else interpolate again
    contour_filt_int_first = contour_filt_first;
    contour_filt_int_last = contour_filt_last;

    % while min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.9
    while ((max(sqrt(diff(real(contour_filt_int_first)).^2+diff(imag(contour_filt_int_first)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int_first)).^2+diff(imag(contour_filt_int_first)).^2))<0.45))
        contour_temp = contour_filt_int_first(:)';
        contour_filt_int_first = Fct_ContourInterpSpline(contour_temp,1,1);
    end

    while ((max(sqrt(diff(real(contour_filt_int_last)).^2+diff(imag(contour_filt_int_last)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int_last)).^2+diff(imag(contour_filt_int_last)).^2))<0.45))
        contour_temp = contour_filt_int_last(:)';
        contour_filt_int_last = Fct_ContourInterpSpline(contour_temp,1,1);
    end

    % Extract the filtered contour points coordinates
    x_filt_int_first = real(contour_filt_int_first);
    y_filt_int_first = abs(imag(contour_filt_int_first));
    x_filt_int_last = real(contour_filt_int_last);
    y_filt_int_last = abs(imag(contour_filt_int_last));

    figure(3)
    subplot(221); imagesc(RawIm_first); colormap gray; title('First image : raw'); hold on; plot(x_filt_int_first,y_filt_int_first,'r-'); hold off; axis off
    subplot(222); imagesc(BinIm_first); title('First image : binarized'); axis off; 
    subplot(223); imagesc(RawIm_last); colormap gray; title('Last image : raw'); hold on; plot(x_filt_int_last,y_filt_int_last,'r-'); hold off; axis off
    subplot(224); imagesc(BinIm_last); title('Last image : binarized'); axis off; 
    pause(3)

    button = questdlg('Threshold OK ?','BinThreshold');
        switch button
            case 'Yes'
                BinThreshold = thresholdtemp;
                threshold_test = false;
            case 'No'
                Para = inputdlg({'New coefficient : '}, 'coeff', 1.0, {'1.0'});
                coeff = str2double(Para{1});
            case 'Cancel'
                break
        end     
end
close(3)

%% BINARISATION

disp('### BINARISATION ###')

tic
for i = nbIm_start:nbIm 
    if HIT
        disp(['Binarisation ' Name_AdjustIm(i,:) '...']);
        AdjustIm = im2double(imread([AdjustDir '\' Name_AdjustIm(i,:)])); % Reading the image to process
        
        BinIm = imbinarize(AdjustIm,BinThreshold); % Binarize the full image with the computed threshold
        BinIm = iminv(BinIm); % Returns the negative of the image
        BinIm = BinIm.*Mask_viewport; % Sets pixel outside the viewport mask to black
    
        % Erode and Dilate the image to fill the gaps and filter the
        % binarisation noise
    
        % BinIm = imerode(BinIm,se); % Erosion and dilatation with mask se to remove binarisation noise
        % BinIm = imdilate(BinIm,se);
        BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure in the image
        BinIm = medfilt2(BinIm, [smfilt smfilt]); % Median filter to remove the last binarisation noise
    
        % BinIm = imdilate(BinIm,se);
        % BinIm = imdilate(BinIm,se);
        % BinIm = imdilate(BinIm,se);
        % BinIm = imdilate(BinIm,se);
        % BinIm = imdilate(BinIm,se);
        % BinIm = imerode(BinIm,se);
        % BinIm = imerode(BinIm,se);
        % BinIm = imerode(BinIm,se);
        % BinIm = imerode(BinIm,se);
        % BinIm = imerode(BinIm,se);
    
        BinIm = imfill(BinIm,'holes'); % Fill the detection gaps in the image
        BinIm = medfilt2(BinIm, [5 5]); % Median filter to remove the last binarisation noise
        BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure on the image (normally the flame)
    
        % figure
        % imshowpair(Im,BinIm,'montage')
        
        imwrite(BinIm,strcat([BinDir '\Binary_'],List_CorrecIm(i).name(end-8:end-4),'.TIF'),'TIF'); % Save the image and process the following 

        % Check if the detected flame exit the HIT zone
        if max(BinIm.*(1-Mask_HIT),[],'all')==1
            HIT = false; % If yes don't process aditionnal images
            disp('Flame exited the HIT zone')         
        end
    end
end

disp('### BINARISATION DONE ###')
toc
%% CONTOUR DETECTION

disp("### CONTOUR DETECTION ###")

List_BinIm = dir([BinDir '\*.tif']); % Gather all binary images from the directory
Name_BinIm = sortrows(char(List_BinIm.name)); % Name of all binary images
nbIm = size(List_BinIm,1); % Number of images to process

tic
for i = 1:nbIm
    
    disp(['Contour detection ' Name_BinIm(i,:) '...']);
    BinIm = im2double(imread([BinDir '\' Name_BinIm(i,:)])); % Reading the binary image
    
    fig_cont = figure('Visible','off');
    cont = imcontour(BinIm,1); % Contour detection
    close(fig_cont)

    x_cont = cont(1,2:end); % Contour coordinates
    y_cont = cont(2,2:end);

    % Check if contour is closed
    if x_cont(end)==x_cont(1) && y_cont(end) == y_cont(1)
        if x_cont(1)== x_cont(end) % If contour extremities are aligned
            x_cont = [x_cont x_cont(1)]; % closing manualy the contour
            y_cont = [y_cont y_cont(1)];
        else
            slope = (y_cont(end)-y_cont(1))/(x_cont(end)-x_cont(1)); % Closing the contour using linear segment
            y_intersect = y_cont(1)-slope*x_cont(1);
            x_line = x_cont(1)+0.5:0.5:x_cont(end)-0.5;
            y_line = slope.*x_line + y_intersect;
            x_cont = [x_cont fliplr(x_line) x_cont(1)];
            y_cont = [y_cont fliplr(y_line) y_cont(1)];
        end
    end

    % Reconstructing the contour geometry
    contour = x_cont -1i * y_cont; 

    % Contour filtering
    contour_filt = Fct_Contour_Filter(contour, ScaleFilter, Magn);

    % Contour interpolation
    % Compute the distance between two successive points 
    % If interpolation is good the distance should be around 1
    % Else interpolate again
    contour_filt_int = contour_filt;
    % while min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.9
    while ((max(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.45))
        contour_temp = contour_filt_int(:)';
        contour_filt_int = Fct_ContourInterpSpline(contour_temp,1,1);
    end

    x_filt_int = real(contour_filt_int);
    y_filt_int = abs(imag(contour_filt_int));

    % Plot the contour on the corresponding raw image for validation
    RawIm = imread([WorkDir '\' Name_RawIm(i,:)]);
    contour_plot = figure('Visible','off');
    imagesc(RawIm); axis equal; title(['Image n°' num2str(i)]); colormap gray; hold on; plot(x_filt_int,y_filt_int,'r-'); plot(cont_HIT(1,2:end),cont_HIT(2,2:end),'r--',LineWidth=2); hold off;
    saveas(contour_plot,[ContourDir '\Image_' Name_RawIm(i,end-8:end-4) '.tif']);
    close(contour_plot)

    save([ContourDir '\Contour_' Name_RawIm(i,end-8:end-4)], 'contour','contour_filt','contour_filt_int','ScaleFilter','Magn','nbIm_start');
    
end
disp('### CONTOUR DETECTION DONE ###')
toc

