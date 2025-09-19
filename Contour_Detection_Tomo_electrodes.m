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
Original_WorkDir = WorkDir; % Keep in memory the original working directory

TruncDir = [WorkDir '\1_Truncation']; % Directory for truncated images
mkdir(TruncDir)
BinDir = [WorkDir '\2_Binary']; % Directory for binary images
mkdir(BinDir)
ContourDir = [WorkDir '\3_Contour']; % Directory for contour data
mkdir(ContourDir)


%% PARAMETERS

load([WorkDir '\Parameters']) % Load processing parameters

HIT_size = 40; % HIT zone size [mm]
HIT_size_pxl = floor(HIT_size/Magn); % HIT zone size [pxl]

Trunc = true; % true = Truncate all images, true by default, put false if images have already been truncated
nbIm_start = 20; % First image to be processed for contour detection

HIT = true; % True by default, control if a flame is within the HIT zone

% Erosion parameters

smfilt = 7; % Size of the median filter
erosion_size = 2;
se = strel('disk',erosion_size,6); % Erosion structure

ScaleFilter = 0.4; % Scale of filtering the pixelisation noise (mm)

%% BACKGROUND IMAGE AND TRUNCATION

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


if Trunc
    disp('### TRUNCATION ###')
    tic 
    for i = 1:nbIm
        disp(['Truncation ' Name_RawIm(i,:) '...']);
        RawIm = imread([WorkDir '\' Name_RawIm(i,:)]);
        TruncIm = (imsubtract(RawIm,Background)); % Substracting the background 
        imwrite(TruncIm,[TruncDir '\Trunc_' num2str(i,'%04.0f') '.tif']); % Saving the truncated image
    end
    disp('### TRUNCATION DONE ###')
    toc
end

%% BACKGROUND IMAGE AND INTENSITY PROFILE

WorkDir = TruncDir; % Switching working directory
List_TruncIm = dir([WorkDir '\*.tif']); % Gather all truncated images in the directory
BACKGR = exist([WorkDir '\mean_background.tif'], 'file'); % check if a mean background image already exists

laserProfile = 0; % Initialising variables
laserProfile_filt = 0;

% Creating a background image as the mean of the first 10 pictures
if BACKGR==0
    meanbackIm = zeros(sizey,sizex); % Empty image
    for i=1:10 % Read the remaining images.
        Im = im2double(imread([WorkDir '\' List_TruncIm(i).name]));
        meanbackIm = meanbackIm + Im;
        delete([WorkDir '\' List_TruncIm(i).name]) % Deleting used image
        % Extracting the laser intensity profile
        laserProfile = laserProfile + mean(Im,2);
        laserProfile_filt = laserProfile_filt + filtfilt(ones(1,100),100,mean(Im,2));
    end
    meanbackIm = meanbackIm / 10; 
    imwrite(meanbackIm,[WorkDir '\mean_background.tif']); % Create background.tif based on the first picture
    laserProfile = laserProfile./10;
    laserProfile_filt = laserProfile_filt./10;
    save([Original_WorkDir '\Laser_Profile.mat'],'laserProfile', 'laserProfile_filt');
    clear meanbackIm Im
end

% Adjusting which images are processed
List_TruncIm = dir([WorkDir '\*.tif']);  % Gather all raw images in the directory
List_TruncIm = List_TruncIm(1:end-1); % removing mean_background.tif
nbIm = size(List_TruncIm,1); % Number of images to process
Name_TruncIm = sortrows(char(List_TruncIm.name)); % Name of the images
meanbackIm = im2double(imread([WorkDir '\mean_background.tif'])); % Reading the mean background image
load([Original_WorkDir '\Laser_Profile.mat'])

% % Extracting the laser intensity profile from the background image
% laserProfile = mean(meanbackIm,2); % Average lines intensity for the background image
% laserProfile_filt = filtfilt(ones(1,100),100,laserProfile); % Filtering the intensity profile

% Plot for visual check 
figure(1)
hold on
plot(laserProfile)
plot(laserProfile_filt,'-r')
hold off
legend('Laser profile','Filtered laser profile')
saveas(gcf, [Original_WorkDir '\Laser_Profile.fig']);

%% MASKS

load electrodes_mask_L.mat % Loading default masks
load electrodes_mask_R.mat
load viewport_mask.mat 
figure,imshow(imread([Original_WorkDir '\' Name_RawIm(end,:)])) % Show the last raw image 
h1 = drawpolygon(gca, 'Position', electrodes_posL); % Draw the mask for the electrodes
wait(h1);
disp('#    ELECTRODES LEFT MASK OK     #')
figure,imshow(imread([Original_WorkDir '\' Name_RawIm(end,:)]))
h2 = drawpolygon(gca, 'Position', electrodes_posR); % Draw the mask for the electrodes
wait(h2);
disp('#    ELECTRODES RIGHT MASK OK     #')
figure,imshow(imread([Original_WorkDir '\' Name_RawIm(end,:)])) 
h3 = drawpolygon(gca, 'Position', viewport_pos); % Draw the mask for the viewport
wait(h3); 
disp('#    VIEWPORT MASK OK     #')

Mask_electrodeL = createMask(h1); % Create a binary mask for the electrodes
Mask_electrodeR = createMask(h2);
Mask_viewport = createMask(h3); % Create a binary mask for the viewport

electrodes_posL = h1.Position; % Extracting position
electrodes_posR = h2.Position;
viewport_pos = h3.Position;

save('electrodes_mask_L.mat','electrodes_posL') % Save the created masks
save('electrodes_mask_R.mat','electrodes_posR')
save('viewport_mask.mat','viewport_pos')

% Finding the center of the 2 electrodes
tip_electrodeL = find(electrodes_posL(:,1)==max(electrodes_posL(:,1)));
tip_electrodeR = find(electrodes_posR(:,1)==min(electrodes_posR(:,1)));
x_electrodes_center = floor(electrodes_posL(tip_electrodeL,1) + (electrodes_posR(tip_electrodeR,1)-electrodes_posL(tip_electrodeL,1))/2);
y_electrodes_center = floor(electrodes_posL(tip_electrodeL,2) + (electrodes_posR(tip_electrodeR,2)-electrodes_posL(tip_electrodeL,2))/2);
x_min_HIT = x_electrodes_center - HIT_size_pxl/2;
y_min_HIT = y_electrodes_center - HIT_size_pxl/2;

button = questdlg('Continuous flame ?');
    switch button
        case 'Yes'
            continuous_flame = true;
        case 'No'
            continuous_flame = false;
    end

close all

fig = figure('Visible','off');
imshow(imread([Original_WorkDir '\' Name_RawIm(end,:)]))
h4 = drawrectangle(gca,'Position',[x_min_HIT,y_min_HIT,HIT_size_pxl,HIT_size_pxl]); % Draw a rectangle corresponding to the HIT zone
Mask_HIT = createMask(h4); % Create a binary mask for the HIT zone

close all

%% BINARY THRESHOLD 

threshold_test = true; % Controls the loop
coeff = 1.0; % Coefficient to adjust the threshold

RawIm_first = imread([Original_WorkDir '\' Name_RawIm(10+nbIm_start,:)]);
RawIm_last = imread([Original_WorkDir '\' Name_RawIm(end,:)]);
TruncIm_first = im2double(imread([WorkDir '\' Name_TruncIm(nbIm_start,:)]));
TruncIm_last = im2double(imread([WorkDir '\' Name_TruncIm(end,:)]));

% Normalisation to take account of the non homogeneity of the laser sheet
for coll = 1:sizey
        TruncIm_first(:,coll) = TruncIm_first(:,coll)./(laserProfile_filt/max(laserProfile_filt)); 
        TruncIm_last(:,coll) = TruncIm_last(:,coll)./(laserProfile_filt/max(laserProfile_filt)); 
end

thresholdtemp = graythresh(TruncIm_first); % temporary threshold used for the loop
fprintf('Binary threshold = %0.03f', thresholdtemp); % Affiche la valeur de seuil trouvée

while threshold_test
    thresholdtemp = thresholdtemp*coeff; % Updating the threshold

    % Full binarisation process for the first and last images
    BinIm_first = imbinarize(TruncIm_first, thresholdtemp);
    BinIm_first = iminv(BinIm_first); 
    BinIm_first = BinIm_first.*Mask_viewport;
    if continuous_flame
    BinIm_first = Fct_Struct_Max(BinIm_first); 
    else
        BinIm_first_label = bwlabel(BinIm_first,4);
        nb_struct = max(BinIm_first_label(:));
        for num_struct = 1:nb_struct
            Pix_in_struct = find(BinIm_first_label==num_struct);
            if length(Pix_in_struct)<200
                BinIm_first(Pix_in_struct)=0;
            end
        end
    end
    BinIm_first = BinIm_first.*iminv(Mask_electrodeL); 
    BinIm_first = BinIm_first.*iminv(Mask_electrodeR);
    BinIm_first = imerode(BinIm_first,se); 
    BinIm_first = imdilate(BinIm_first,se);
    %BinIm_first = Fct_Struct_Max(BinIm_first); 
    BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]);
    BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = imdilate(BinIm_first,se);
    BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imerode(BinIm_first,se);
    BinIm_first = imfill(BinIm_first,'holes'); 
    BinIm_first = Fct_Struct_Max(BinIm_first); 
    BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]); 

    BinIm_last = imbinarize(TruncIm_last, thresholdtemp);
    BinIm_last = iminv(BinIm_last); 
    BinIm_last = BinIm_last.*Mask_viewport;
    if continuous_flame
    BinIm_last = Fct_Struct_Max(BinIm_last); 
    else
        BinIm_last_label = bwlabel(BinIm_last,4);
        nb_struct = max(BinIm_last_label(:));
        for num_struct = 1:nb_struct
            Pix_in_struct = find(BinIm_last_label==num_struct);
            if length(Pix_in_struct)<200
                BinIm_last(Pix_in_struct)=0;
            end
        end
    end
    BinIm_last = BinIm_last.*iminv(Mask_electrodeL); 
    BinIm_last = BinIm_last.*iminv(Mask_electrodeR);
    BinIm_last = imerode(BinIm_last,se); 
    BinIm_last = imdilate(BinIm_last,se);
    %BinIm_last = Fct_Struct_Max(BinIm_last); 
    BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]);
    BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = imdilate(BinIm_last,se);
    BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imerode(BinIm_last,se);
    BinIm_last = imfill(BinIm_last,'holes'); 
    BinIm_last = Fct_Struct_Max(BinIm_last); 
    BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]); 

    figure(3)
    subplot(221); imagesc(RawIm_first); colormap gray; title('First image : raw'); axis off;
    subplot(222); imagesc(BinIm_first); title('First image : binarized'); axis off; 
    subplot(223); imagesc(RawIm_last); title('Last image : raw'); axis off; 
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
        disp(['Binarisation ' Name_TruncIm(i,:) '...']);
        Im = im2double(imread([WorkDir '\' Name_TruncIm(i,:)])); % Reading the image to process
        for coll = 1:sizey
            Im(:,coll) = Im(:,coll)./(laserProfile_filt/max(laserProfile_filt)); % Normalisation to take account of the non-homogeneity of the laser sheet
        end
        
        Im2 = Im./Mask_viewport; % Applying the viewport mask on a temporary image
    %     Index_Visu = isinf(Im2); % Find each inf pixel in Im2 ie. outside the mask
    %     BinThreshold = graythresh(Im2(~Index_Visu)); % Compute the gray threshold for the binarization using Otsu method
        BinIm = imbinarize(Im,BinThreshold); % Binarize the full image with the computed threshold
        BinIm = iminv(BinIm); % Returns the negative of the image
        BinIm = BinIm.*Mask_viewport; % Sets pixel outside the viewport mask to black
        if continuous_flame
            BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure in the image
        else
            BinIm_label = bwlabel(BinIm,4);
            nb_struct = max(BinIm_label(:));
            for num_struct = 1:nb_struct
                Pix_in_struct = find(BinIm_label==num_struct);
                if length(Pix_in_struct)<200
                    BinIm(Pix_in_struct)=0;
                end
            end
        end
        BinIm = BinIm.*iminv(Mask_electrodeL); % Removing electrodes using the masks
        BinIm = BinIm.*iminv(Mask_electrodeR);
    
        % Erode and Dilate the image to fill the gaps and filter the
        % binarisation noise
    
        BinIm = imerode(BinIm,se); % Erosion and dilatation with mask se to remove binarisation noise
        BinIm = imdilate(BinIm,se);
        % BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure in the image
        BinIm = medfilt2(BinIm, [smfilt smfilt]); % Median filter to remove the last binarisation noise
    
        BinIm = imdilate(BinIm,se);
        BinIm = imdilate(BinIm,se);
        BinIm = imdilate(BinIm,se);
        BinIm = imdilate(BinIm,se);
        BinIm = imdilate(BinIm,se);
        BinIm = imerode(BinIm,se);
        BinIm = imerode(BinIm,se);
        BinIm = imerode(BinIm,se);
        BinIm = imerode(BinIm,se);
        BinIm = imerode(BinIm,se);
    
        BinIm = imfill(BinIm,'holes'); % Fill the detection gaps in the image
        BinIm = medfilt2(BinIm, [5 5]); % Median filter to remove the last binarisation noise
        BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure on the image (normally the flame)
    
        % figure
        % imshowpair(Im,BinIm,'montage')

        % Check if the detected flame exit the HIT zone
        if max(BinIm.*(1-Mask_HIT),[],'all')==1
            HIT = false; % If yes don't save and don't process aditionnal images
            disp('Flame exited the HIT zone')
        else
            imwrite(BinIm,strcat([BinDir '\Binary_'],List_TruncIm(i).name(end-6:end-4),'.TIF'),'TIF'); % If no save the image and process the following 
        end
    end
end

disp('### BINARISATION DONE ###')
toc
%% CONTOUR DETECTION

disp("### CONTOUR DETECTION ###")

WorkDir = BinDir;
List_BinIm = dir([WorkDir '\*.tif']); % Gather all binary images from the directory
Name_BinIm = sortrows(char(List_BinIm.name)); % Name of all binary images
nbIm = size(List_BinIm,1); % Number of images to process

tic
for i = 1:nbIm
    
    disp(['Contour detection ' Name_BinIm(i,:) '...']);
    BinIm = im2double(imread([WorkDir '\' Name_BinIm(i,:)])); % Reading the binary image
    
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

    % Remove points in electrode mask from the contour
    % for contour_point = 1:size(x_cont,2)
    %     % Left electrode
    %     [in, on] = inpolygon(x_cont(contour_point),y_cont(contour_point),electrodes_posL(:,1),abs(electrodes_posL(:,2)));
    %     if in || on 
    %         x_cont(contour_point) = NaN; % replacing values by NaN to keep the array size
    %         y_cont(contour_point) = NaN;
    %     end
    %     % Right electrode
    %     [in, on] = inpolygon(x_cont(contour_point),y_cont(contour_point),electrodes_posR(:,1),abs(electrodes_posR(:,2)));
    %     if in || on 
    %         x_cont(contour_point) = NaN;
    %         y_cont(contour_point) = NaN;
    %     end
    % end
    x_cont = x_cont(1,~isnan(x_cont(1,:)));
    y_cont = y_cont(1,~isnan(y_cont(1,:)));

    contour = x_cont -1i * y_cont; % Reconstructing the contour geometry

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
    % RawIm = imread([Original_WorkDir '\' Name_RawIm(i+nbIm_start+9,:)]);
    RawIm = imread([WorkDir '\' Name_BinIm(i,:)]);
    contour_plot = figure('Visible','off');
    imagesc(RawIm); axis equal; title(['Image n°' num2str(i)]); colormap gray; hold on; plot(x_filt_int,y_filt_int,'r-');
    saveas(contour_plot,[ContourDir '\Image_' Name_BinIm(i,end-6:end-4) '.tif']);
    close(contour_plot)

    save([ContourDir '\Contour_' Name_BinIm(i,end-6:end-4)], 'contour','contour_filt','contour_filt_int','ScaleFilter','Magn','nbIm_start','F');
    
end
disp('### CONTOUR DETECTION DONE ###')
toc

