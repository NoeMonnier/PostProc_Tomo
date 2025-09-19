%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PROPAGATION SPEED TOMO                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% WORKING DIRECTORY


addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\'); % Working directory 
Original_WorkDir = WorkDir; % Keep in memory the original working directory

ContourDir = [WorkDir '\3_Contour']; % Directory for contour data
SpeedDir = [WorkDir '\4_Propagation_Speed']; % Directory for propagation speed data
mkdir(SpeedDir)

% load electrodes_mask_L.mat % Load the mask coordinates 
% load electrodes_mask_R.mat

%% PARAMETERS

load([WorkDir '\Parameters']) % Load camera parameters
R_min = 5 ; % Smallest processed radius [mm]
dt = 1/F*1000; % Timestep between 2 images [ms] 

grad_size = 1; % 5 by default

Speed_color=[0 0.5 0]; % Color array for local speed plot

%% CONTOUR GEOMETRY

List_Contours = dir([ContourDir '\*.mat']); % Gather all the contour data
nbContours = size(List_Contours,1); % Number of contours
NumContours = zeros(nbContours,1); % Store the numbers of processed contours

% Pre allocating arrays to stock computed data
A = zeros(nbContours,1); % Flame contour area
P = zeros(nbContours,1); % Flame contour perimeter
X_center = zeros(nbContours,1); % Flame contour center coordinates
Y_center = zeros(nbContours,1);
Ra = zeros(nbContours,1); % Area based radii
Rp = zeros(nbContours,1); % Perimeter based radii
Wrink_mean = zeros(nbContours,1); % Flame wrinkeling
twoD_sphericity = zeros(nbContours,1); % 2D-Sphericity
roundness = zeros(nbContours,1); % Roundness

for i = 1:nbContours
    
    load([ContourDir '\' List_Contours(i).name]) % Load the contour 
    NumContour = List_Contours(i).name(end-8:end-4); % Find the number of the contour
    NumContours(i) = str2double(NumContour); % Store the contour numbers as a double array

    XContour = real(contour_filt_int); % Gather de coordinates of every points in the contour
    YContour = abs(imag(contour_filt_int));

    [geom, ~,~] = polygeom(XContour,YContour); % Compute the geometric parameters of the contour
    A(i) = geom(1)*(Magn^2); % Area of the flame in mm²
    P(i) = geom(4)*Magn; % Perimeter of the flame in mm
    X_center(i) = geom(2)*Magn; % X coordinate of the flame center
    Y_center(i) = geom(3)*Magn; % Y coordinate of the flame center

    % Radii are computed based on a spherical flame assumption (false)
    Ra(i) = sqrt(A(i)/pi()); % Area based radius in mm
    Rp(i) = P(i)/(2*pi()); % Perimeter based radius in mm
    Wrink_mean(i) = (Rp(i)^2)/(Ra(i)^2); % Wrinkling of the flame 

    % Spherical approching parameters are computed
    twoD_sphericity(i) = (2*pi()*sqrt(A(i)/pi()))/P(i); % 2D sphericity (closer to one for a shape approaching a 2D-sphere) 
    roundness(i) = (4*pi()*A(i))/(P(i)^2); % Roundness (close to 1 = close to a circle)

end

% Reconstructing time array in ms
time(:,1) = (NumContours./20).*dt; 

%% GLOBAL SPEED

% for i = 2:nbContours
% 
%     Vt_A(i) = (A(i)-A(i-1))/P(i)/dt; % Global propagation speed in m/s
%     Stretch_global(i) = (A(i)-A(i-1))/A(i)/dt; % Global stretch rate in 1/ms
% 
% end

A = A*1e-6; % (mm²) ==> (m²)
P = P*1e-3; % (mm) ==> (m)
dAdt = grad(A,dt/1000,grad_size);
Vt_A = dAdt./P; % global propagation speed based on area growth [m/s]
dR_Adt = grad(Ra,dt,grad_size); % Turbulent propagation speed based on area radius [m/s]
Stretch_global = dAdt./A; % Global stretch rate of the flame [1/s]

%% LOCAL SPEED ALONG THE FLAME CONTOUR

i_start = 12; % Starting index for the discretisation
nb_points = 150; % Number of discrete points in the contour

Contour_Rmin = find(Ra>R_min,1); % Find the 1st contour with radius greater than the minimal threshold
disp(['1st processed contour : ' num2str(Contour_Rmin)])
nbContours_processed = nbContours-Contour_Rmin; 
disp(['Number of processed contour : ' num2str(nbContours_processed)])

matrix_index = 1; % Initialisation of the index for the matrix storing solution
tic
for i = Contour_Rmin:nbContours-1

    Contour1 = load([ContourDir '\' List_Contours(i).name]); % Load both i and i+1 contours
    Contour2 = load([ContourDir '\' List_Contours(i+1).name]);

    NumContour1 = List_Contours(i).name(end-6:end-4); % Numbers of the studied contours
    NumContour2 = List_Contours(i+1).name(end-6:end-4);

    disp(['Local speeds computation for contour ' NumContour1 ])

    XContour1 = real(Contour1.contour_filt_int); 
    YContour1 = abs(abs(imag(Contour1.contour_filt_int))-sizey); 
    XContour2 = real(Contour2.contour_filt_int);
    YContour2 = abs(abs(imag(Contour2.contour_filt_int))-sizey);

    step_contour = floor((length(XContour1)-12-i_start)/nb_points); % Step size between 2 descreate points
    if length(XContour1)<=step_contour*nb_points % If step size is too big reduce step size
      nb_points=length(XContour1);
      pas_contour=1;
    end

    solution_index = 1; % initialisation of the solution index for local speed computation
    nbContour_point = length(XContour1)-12-i_start; % Number of contour points to be computed 
    % Mermory allocation for solutions arrays
    Xfinal_contour1 = zeros(1,nbContour_point);
    Yfinal_contour1 = zeros(1,nbContour_point);
    Xfinal_contour2 = zeros(1,nbContour_point);
    Yfinal_contour2 = zeros(1,nbContour_point);
    distance_x = zeros(1,nbContour_point);
    distance_y = zeros(1,nbContour_point);
    Speed_U_pix = zeros(1,nbContour_point);
    Speed_V_pix = zeros(1,nbContour_point);
    Speed_local = zeros(1,nbContour_point);

    % Loop on all descreate point to compute local speed
    for contour_point = i_start:step_contour:length(XContour1)-12

        tangX = XContour1(contour_point+3)-XContour1(contour_point-3); % Horizontal tangeante 
        tangY = YContour1(contour_point+3)-YContour1(contour_point-3); % Vertical tangeante
        % Tangeante equation at contour_point coordinates
        b_tang = YContour1(contour_point)-tangY/tangX*XContour1(contour_point);
        Xcoordinate = XContour1(contour_point)-4:0.2:XContour1(contour_point)+4;
        y_tang = tangY/tangX*Xcoordinate + b_tang;
        % Normal equation at contour_point coordinates
        b_norm = YContour1(contour_point)+tangX/tangY*XContour1(contour_point);
        y_norm = -tangX/tangY*Xcoordinate+b_norm;

        % We search a point on the 2nd contour that is closed to the normal
        interval_gap = 10;
        % list of all points close to the studied point 
        list_points_contour2 = find(abs(XContour2-XContour1(contour_point))<=interval_gap & abs(YContour2-YContour1(contour_point))<=interval_gap);
        % If list is empty increase the interval gap and recompute the list
        while isempty(list_points_contour2)
            interval_gap = interval_gap+2;
            list_points_contour2 = find(abs(XContour2-XContour1(contour_point))<=interval_gap & abs(YContour2-YContour1(contour_point))<=interval_gap);
        end

        if length(list_points_contour2)<=1
            contour2_point=list_points_contour2; % If only one point is within the interval we select it        
        else
            nb_index=floor(length(list_points_contour2)/2); % If more than 1 point is in the interval we take the middle point 
            contour2_point=list_points_contour2(nb_index); % index of the point on the 2nd contour
        end
        
        % If the selected point is close from the first or the last point
        % of the contour, we reduce the research interval
        zone_search = 10;
        while contour2_point-zone_search<1 || contour2_point+zone_search>length(XContour2)
            zone_search=zone_search-1;
        end

        % Interpolation of the contour and the normal as polynomes to solve
        % the geometric equation contour = normale and find the
        % intersection point

        if contour2_point>10

            % We only use the point in the search zone for the
            % interpolation
            XContour2_zone=XContour2(contour2_point-zone_search:1:contour2_point+zone_search);
            Ycontour2_zone=YContour2(contour2_point-zone_search:1:contour2_point+zone_search);
    
            coeff_poly_contour=polyfit(XContour2_zone,Ycontour2_zone,2); % Polynomial coefficients for the contour (2nd order polynome)
            coeff_poly_norm=polyfit(Xcoordinate,y_norm,1); % Polynomial coefficients for the normal (1st order polynome)

            % Refining the interpolation zone
            if XContour2(contour2_point-zone_search)>=XContour2(contour2_point+zone_search)
                XContour2_zone=XContour2(contour2_point-zone_search):-0.1:XContour2(contour2_point+zone_search); % Research zone on contour 2
            else
                XContour2_zone=XContour2(contour2_point-zone_search):0.1:XContour2(contour2_point+zone_search);
            end
            XContour1_zone = XContour1(contour_point)-5:0.1:XContour1(contour_point)+5; % Research zone on contour 1
        else

            XContour2_zone=XContour2(1:1:contour2_point+zone_search);
            YContour2_zone=YContour2(1:1:contour2_point+zone_search);
            coeff_poly_contour=polyfit(XContour2_zone,YContour2_zone,2); % Polynomial coefficients for the contour (2nd order polynome)
            coeff_poly_norm=polyfit(Xcoordinate,y_norm,1); % Polynomial coefficients for the normal (1st order polynome)
    
            if XContour2(1)>=XContour2(contour2_point+zone_search)
                XContour2_zone=XContour2(1):-0.1:XContour2(contour2_point+zone_search);
            else
                XContour2_zone=XContour2(1):0.1:XContour2(contour2_point+zone_search);    
            end
                XContour1_zone=XContour1(contour_point)-5:0.1:XContour1(contour_point)+5;
        end

        poly_contour2 = polyval(coeff_poly_contour,XContour2_zone); % Polynomial interpolation of contour 2
        poly_norm = coeff_poly_norm(1).*XContour1_zone + coeff_poly_norm(1); % Polynomial interpolation of contour normal

        % Solving the system : ax+b = c1x^2 + c2x + c3
        % Solutions represents the intersections between the normal and
        % contour2
        syms x
        solutions_x=eval(solve(coeff_poly_contour(1).*x.^2+coeff_poly_contour(2).*x-coeff_poly_norm(1).*x-coeff_poly_norm(2)+coeff_poly_contour(3), x));
        % Complex solutions are discarded
        for k=1:length(solutions_x)
            if imag(solutions_x(k))>0
                solutions_x(k)=NaN;
            end
        end
        
        solution_criterion = 1; % Criterion to discriminate between 2 real solutions
        finalX_index = find((abs(solutions_x-XContour2(contour2_point)))<solution_criterion); % Find the index of x corresponding to the solution
        % If no solutions is found, increase the criterion 
        while isempty(finalX_index) 
            solution_criterion=solution_criterion+1;
            finalX_index=find(abs(solutions_x-XContour2(contour2_point))<solution_criterion); 
        end

        % If 2 solutions are found discriminate using the minimum variation 
        if length(finalX_index)==2
            [~, finalX_index]=min(abs(solutions_x(finalX_index)-XContour2(contour2_point)));
        end
        
        % Compute the final coordinate of the selected point on the 2
        % contours
        Xfinal_contour1(solution_index)=XContour1(contour_point);
        Yfinal_contour1(solution_index)=YContour1(contour_point);
        Xfinal_contour2(solution_index)=solutions_x(finalX_index);
        Yfinal_contour2(solution_index)=coeff_poly_contour(1).*Xfinal_contour2(solution_index).^2+coeff_poly_contour(2).*Xfinal_contour2(solution_index)+coeff_poly_contour(3);

        % Distance and speed are computed and stored for plots
        distance_x(solution_index)=Xfinal_contour2(solution_index)-Xfinal_contour1(solution_index);
        distance_y(solution_index)=Yfinal_contour2(solution_index)-Yfinal_contour1(solution_index);
        Speed_U_pix(solution_index)=distance_x(solution_index)/dt;
        Speed_V_pix(solution_index)=distance_y(solution_index)/dt;

        % If the points are too far away due to flame curvature being too
        % high, the solution is discarded
        if abs(distance_x(1,solution_index))>=25 || abs(distance_y(1,solution_index))>=25
            Xfinal_contour2(solution_index)=NaN;
            Yfinal_contour2(solution_index)=NaN;
            distance_x(solution_index)=NaN;
            distance_y(solution_index)=NaN;
            Speed_U_pix(solution_index)=NaN;
            Speed_V_pix(solution_index)=NaN;
        end
        
        % Computation of the local speed (norm of speed vector)
        Speed_local(solution_index) = sqrt((Speed_U_pix(solution_index))^2+(Speed_V_pix(solution_index))^2)*Magn*1e-3;
        % If the norm is complex the solution is discarded
        if abs(imag(Speed_local(solution_index)))>0
           Speed_local(solution_index)=NaN;
        end

        % check if the computed points are in the electrodes masks, if yes
        % the data will be ignore to avoid introducing local error 
        % [in, on] = inpolygon(XContour1(contour_point),YContour1(contour_point),electrodes_posL(:,1),abs(electrodes_posL(:,2)-sizey));
        % if in || on 
        %     Speed_U_pix(solution_index) = NaN;
        %     Speed_V_pix(solution_index) = NaN;
        %     Speed_local(solution_index) = NaN;
        % end
        % [in, on] = inpolygon(XContour1(contour_point),YContour1(contour_point),electrodes_posR(:,1),abs(electrodes_posR(:,2)-sizey));
        % if in || on 
        %     Speed_U_pix(solution_index) = NaN;
        %     Speed_V_pix(solution_index) = NaN;
        %     Speed_local(solution_index) = NaN;
        % end

        solution_index = solution_index + 1; % Iterate to the next contour point
        % Clearing the array for next points to avoid data overflow
%         clear coeff_poly_contour coeff_poly_norm XContour2_zone YContour2_zone solutions_x poly_contour2 
%         clear list_points_contour2 nb_index contour2_point interval_gap solution_criterion
%         clear x
    end

    LocalSpeed_matrix(matrix_index,1:length(Speed_local))=NaN;
    LocalSpeed_matrix(matrix_index,1:length(Speed_local))=Speed_local(:);
    clear Speed_local
    matrix_index=matrix_index+1; % Iterate to the next contour

    % Plot for the local speed evolution
    LocalSpeed_fig=figure('Visible','off');
    % Contours are centered and plotted
    hold on
    plot(XContour1-X_center(i)/Magn,YContour1-Y_center(i)/Magn,'r-')
    plot(XContour2-X_center(i+1)/Magn,YContour2-Y_center(i+1)/Magn,'b-')
    %plot(electrodes_posL(:,1)-X_center(i)/Magn,abs(electrodes_posL(:,2)-sizey)-Y_center(i)/Magn)
    %plot(electrodes_posR(:,1)-X_center(i)/Magn,abs(electrodes_posR(:,2)-sizey)-Y_center(i)/Magn)
    quiver(Xfinal_contour1(1:2:end)-X_center(i)/Magn, Yfinal_contour1(1:2:end)-Y_center(i)/Magn,Speed_U_pix(1:2:end),Speed_V_pix(1:2:end),1,'Color',Speed_color,'LineWidth',0.8,'Marker','.')
    hold off
    title(['Flame n°' NumContour1 ' evolution'],'FontSize', 12, 'FontName','Times New Roman')
    xlabel('X pixels', 'FontSize',14,'FontName','Times New Roman')
    ylabel('Y pixels', 'FontSize',14,'FontName','Times New Roman')  
    legend([ 'Contour n°' NumContour1], ['Contour n°' NumContour2],'Local speeds', 'Location', 'SouthWest')
    axis([-400 400 -400 400])
    saveas(LocalSpeed_fig,[SpeedDir '\Local_speed' NumContour1 '.tif'],'tif') 
    %saveas(LocalSpeed_fig,[SpeedDir '\Local_speed' NumContour '.fig'],'fig')
    close(LocalSpeed_fig)

    % Clear computing data to avoid data overflow
    clear Xfinal_contour1 Yfinal_contour1 Xfinal_contour2 Yfinal_contour2 Speed_U_pix Speed_V_pix finalX_index solution_index
    clear Contour1 Contour2 XContour1 XContour2 YContour1 YContour2 
end
toc

%% LOCAL SPEED PDF

for i = 1:matrix_index-1

    

end

%% SAVING DATA & RESULTS PLOTS

disp('### Saving radii and speeds ###')
file = [WorkDir '\Results_Rp_Ra_Vt.dat'];
fid = fopen(file, 'w');
fprintf(fid, 'Time (ms)\tPerimeter (m)\tArea (m^2)\tRp (mm)\tRa (mm)\tVt_A (m/s)\tdR_Adt (m/s)\tGlobal Stretch (1/s)\tMean Wrinkling (-)\tX_center (mm)\tY_center (mm)'); 
fprintf(fid,'\n'); % Saut de ligne 
for k = 1:length(time)
    fprintf(fid,'%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f',time(k), P(k), A(k), Rp(k), Ra(k), Vt_A(k), dR_Adt(k), Stretch_global(k), Wrink_mean(k), X_center(k), Y_center(k));
    fprintf(fid, '\n');
end
fclose(fid);
save([WorkDir '\Results_Rp_Ra_Vt'], 'time', 'P', 'A', 'Rp', 'Ra', 'Vt_A','dR_Adt', 'Stretch_global','Wrink_mean', 'X_center', 'Y_center')

% Radii and wrinkling plot
figure(1)
hold on
yyaxis left
plot(time,Ra,'b+')
plot(time,Rp,'g+')
ylabel('Equivalent radius [mm]','Interpreter','latex')
yyaxis right
plot(time,Wrink_mean,'r--')
ylabel('Mean Wrinkling [-]','Interpreter','latex')
hold off
xlabel('Time [ms]','Interpreter','latex')
legend('Ra','Rp','Wrinkling','Interpreter','latex')
xlim([0 time(end)])
box on;
set(gca,'Fontname','Times','Fontsize',16)
saveas(gcf, [WorkDir '\Radii_evolution.fig']) % Save

% Global speed over stretch plot
figure(2)
plot(Stretch_global, Vt_A, 'bo')
xlabel('Stretch rate [s$^{-1}$]','Interpreter','latex')
ylabel('Propagation speed [m/s]','Interpreter','latex')
box on;
set(gca,'Fontname','Times','Fontsize',16)
saveas(gcf, [WorkDir '\Vt_A_Stretch.fig']) % Save

disp('DONE')

% Sphericity plot
figure(3)
hold on
plot(Ra,twoD_sphericity,'r',LineWidth=2)
plot(Ra,roundness,'b',LineWidth=2)
hold off
xlabel('$R_a$ [mm]','Interpreter','latex')
ylabel('[-]')
legend('2D-Sphericity','Roundness','Interpreter','latex')
box on;
set(gca,'Fontname','Times','Fontsize',16)
saveas(gcf, [WorkDir '\Sphericity.fig']) % Save

