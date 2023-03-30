%%
%Task 5a
%%
%
% Main script that reads controller data and laser data

clear all;
close all;

% Co-ordinates of the ref. nodes
REF = [1920 9470;   % 01
       10012 8179;  % 02
       9770 7590;   % 03
       11405 7228;  % 04
       11275 6451;  % 05
       11628 6384.5;% 06
       11438 4948;  % 07
       8140 8274;   % 08
       8392 8486;   % 09
       3280 2750;   % 10
       7250 2085;   % 11
       9990 1620;   % 12
       7485 3225;   % 13
       9505 3893;   % 14
       9602 4278;   % 15
       10412 4150;  % 16
       4090 7920;   % 17
       8010 5290;   % 18
       8255 6099;   % 19
       7733 6151;   % 20
       7490 6136;   % 21
       7061 5420;   % 22
       7634 5342];  % 23

% LINE SEGMENT MODEL (Indeces in the REF-vector)
LINES = [1 8;       % L01
         9 2;       % L02
         2 3;       % L03
         3 4;       % L04
         4 5;       % L05
         5 6;       % L06
         6 7;       % L07
         17 10;     % L08
         10 12;     % L09
         11 13;     % L10
         12 16;     % L11
         16 15;     % L12
         15 14;     % L13
         19 21;     % L14
         22 18;     % L15
         20 23];    % L16
         
% Control inputs (velocity and steering angle)
CONTROL = load('control_joy.txt');

% Laser data
LD_HEAD   = load('laser_header.txt');
LD_ANGLES = load('laser_angles.txt');
LD_MEASUR = load('laser_measurements.txt');
LD_FLAGS  = load('laser_flags.txt');

[no_inputs co] = size(CONTROL);

% Robots initial position
X(1) = CONTROL(1,4);
Y(1) = CONTROL(1,5);
AA(1) = CONTROL(1,6);
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];
sv = .01;%.5;
sa = .2;%2*0.01754;%.2; %least decimal
st = .05/100; % 1 percentage
scan_idx = 1;
fig_path = figure;
fig_env = figure;
ScanPosIndx = [];

% Plot the line model
figure(fig_env); plot_line_segments(REF, LINES, 1);
     
for kk = 2:no_inputs,
    % Check if we should get a position fix, i.e. if the time stamp of the
    % next laser scan is the same as the time stamp of the control input
    % values
    if LD_HEAD(scan_idx,1) == CONTROL(kk-1,1),
        % Mark the position where the position fix is done - and the size
        % of the position fix to be found
        figure(fig_path);
        hold on; plot(X(kk-1), Y(kk-1), 'ro'); hold off;
        hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'ro'); hold off;
        hold on; plot([X(kk-1) CONTROL(kk-1,4)], [Y(kk-1) CONTROL(kk-1,5)], 'r'); hold off;
        
        % > Log pos (index) there scan is taken
        ScanPosIndx = [ScanPosIndx kk-1]; % Used for plotting
        
        % Get the position fix - Use data points that are ok, i.e. with
        % flags = 0
        DATAPOINTS = find(LD_FLAGS(scan_idx,:) == 0);
        angs = LD_ANGLES(scan_idx, DATAPOINTS);
        meas = LD_MEASUR(scan_idx, DATAPOINTS);
        
        % Plot schematic picture of the Snowhite robot
        alfa = 660;
        beta = 0;
        gamma = -90*pi/180;
        figure(fig_env); plot_line_segments(REF, LINES, 1);
        plot_threewheeled_Laser([X(kk-1) Y(kk-1) AA(kk-1)]', 100, 612, 2, CONTROL(kk-1,5), 150, 50, 680, alfa, beta, gamma, angs, meas, 1);
        

        LINEMODEL = [REF(LINES(:,1),1:2) REF(LINES(:,2),1:2)];
 % Task 5a Kalman Filter with simulated position fixes (small) (Exercise 4)
        xv = 10;
        yv = 10;
        asv = 1*pi/180;
        C = [xv^2 0 0;
           0 yv^2 0;
           0 0 asv^2];
        x1 = randn;
        y1 = randn;
        a1 = randn;

        x2 = 0 + xv*x1;
        y2 = 0 + yv*y1;
        a2 = 0 + asv*a1;
        X_1 = CONTROL(kk-1,4) + x2;
        Y_1 = CONTROL(kk-1,5) + y2;
        A_1 = mod(CONTROL(kk-1,6) + a2, 2*pi);
        P1 = P(kk-1,1:9);      
    
        ICP = inv(C + reshape(P1,[3,3]));

        XYZ = [(X_1 -X(kk-1));
            (Y_1 -Y(kk-1));
            (AngDifference(A_1, AA(kk-1)))];
        
        XYZ__ = (reshape(P1,[3,3])*ICP)*XYZ;
        
        
        X(kk-1) = X(kk-1) + XYZ__(1,1);
        Y(kk-1) = Y(kk-1) + XYZ__(2,1);
        AA(kk-1) = mod(AA(kk-1) + XYZ__(3,1), 2*pi);
        
        c_XYZ = inv(inv(C) + inv(reshape(P1,[3,3])));
        
        P(kk-1,1:9) = reshape(c_XYZ,[1,9]);

 

        
        % Next time use the next scan
        scan_idx = mod(scan_idx, max(size(LD_HEAD))) + 1;
    end;
    
    % Mark the estimated (dead reckoning) position
    figure(fig_path);
    hold on; plot(X(kk-1), Y(kk-1), 'b.'); hold off;
    % Mark the true (from the LaserWay system) position
    hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'k.'); hold off;
    
    % Estimate the new position (based on the control inputs) and new
    % uncertainty
    V = CONTROL(kk-1,2);
    A = CONTROL(kk-1,3);
    T = 0.050;
    L = 680;

    X(kk) = X(kk-1) + V*cos(A)*T*cos(AA(kk-1) + (V*sin(A)*T/(2*L)));
    Y(kk) = Y(kk-1) + V*cos(A)*T*sin(AA(kk-1) + (V*sin(A)*T/(2*L)));
    AA(kk) = mod(AA(kk-1) + V*sin(A)*T/L, 2*pi);
    
    % ALSO UPDATE THE UNCERTAINTY OF THE POSITION
    % Task 2 Dead Reckoning - Add you code here


    % Task 2 Dead Reckoning - Add you code here
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];      
    Axya = [1 0 -T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 1  T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 0 1];
    Au = [(T*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T^2*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))  (- T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A) - (T^2*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L))  (V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
        (T*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T^2*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))   ((T^2*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L) - T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A))  (V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
                                                                                        (T*sin(A))/(L)                                                                                (T*V*cos(A))/(L)                                                                                 (V*sin(A))/(L)];
    
    Cu = [sv*abs(V) 0 0;
        0 sa^2 0;
        0 0 st^2];    
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Axya*Cxya_old*Axya' + Au*Cu*Au';

    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];    
end

ERROR = [X' Y' AA'] - CONTROL(:,4:6);
ERROR(:,3) = AngDifference(AA',CONTROL(:,6));
ERROR = abs(ERROR);
figure,
subplot(3,1,1);
plot(ERROR(:,1),'b'); hold;
plot(sqrt(P(:,1)),'r'); % one sigma
plot(ScanPosIndx,sqrt(P(ScanPosIndx,1)),'k.'); % one sigma
title('Error X [mm] and uncertainty [std] (red)');

subplot(3,1,2);
plot(ERROR(:,2),'b'); hold;
plot(sqrt(P(:,5)),'r'); % one sigma
title('Error Y [mm] and uncertainty [std] (red)');

subplot(3,1,3);
plot(ERROR(:,3)*180/pi,'b'); hold;
plot(sqrt(P(:,9))*180/pi,'r'); % one sigma
title('Error A [degree] and uncertainty [std] (red)');
%%
%%
% Task 5b
%%
%
% Main script that reads controller data and laser data

clear all;
close all;

% Co-ordinates of the ref. nodes
REF = [1920 9470;   % 01
       10012 8179;  % 02
       9770 7590;   % 03
       11405 7228;  % 04
       11275 6451;  % 05
       11628 6384.5;% 06
       11438 4948;  % 07
       8140 8274;   % 08
       8392 8486;   % 09
       3280 2750;   % 10
       7250 2085;   % 11
       9990 1620;   % 12
       7485 3225;   % 13
       9505 3893;   % 14
       9602 4278;   % 15
       10412 4150;  % 16
       4090 7920;   % 17
       8010 5290;   % 18
       8255 6099;   % 19
       7733 6151;   % 20
       7490 6136;   % 21
       7061 5420;   % 22
       7634 5342];  % 23

% LINE SEGMENT MODEL (Indeces in the REF-vector)
LINES = [1 8;       % L01
         9 2;       % L02
         2 3;       % L03
         3 4;       % L04
         4 5;       % L05
         5 6;       % L06
         6 7;       % L07
         17 10;     % L08
         10 12;     % L09
         11 13;     % L10
         12 16;     % L11
         16 15;     % L12
         15 14;     % L13
         19 21;     % L14
         22 18;     % L15
         20 23];    % L16
         
% Control inputs (velocity and steering angle)
CONTROL = load('control_joy.txt');

% Laser data
LD_HEAD   = load('laser_header.txt');
LD_ANGLES = load('laser_angles.txt');
LD_MEASUR = load('laser_measurements.txt');
LD_FLAGS  = load('laser_flags.txt');

[no_inputs co] = size(CONTROL);

% Robots initial position
X(1) = CONTROL(1,4);
Y(1) = CONTROL(1,5);
AA(1) = CONTROL(1,6);
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];
sv = .01;%.5;
sa = .2;%2*0.01754;%.2; %least decimal
st = .05/100; % 1 percentage
scan_idx = 1;
fig_path = figure;
fig_env = figure;
ScanPosIndx = [];

% Plot the line model
figure(fig_env); plot_line_segments(REF, LINES, 1);
     
for kk = 2:no_inputs,
    % Check if we should get a position fix, i.e. if the time stamp of the
    % next laser scan is the same as the time stamp of the control input
    % values
    if LD_HEAD(scan_idx,1) == CONTROL(kk-1,1),
        % Mark the position where the position fix is done - and the size
        % of the position fix to be found
        figure(fig_path);
        hold on; plot(X(kk-1), Y(kk-1), 'ro'); hold off;
        hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'ro'); hold off;
        hold on; plot([X(kk-1) CONTROL(kk-1,4)], [Y(kk-1) CONTROL(kk-1,5)], 'r'); hold off;
        
        % > Log pos (index) there scan is taken
        ScanPosIndx = [ScanPosIndx kk-1]; % Used for plotting
        
        % Get the position fix - Use data points that are ok, i.e. with
        % flags = 0
        DATAPOINTS = find(LD_FLAGS(scan_idx,:) == 0);
        angs = LD_ANGLES(scan_idx, DATAPOINTS);
        meas = LD_MEASUR(scan_idx, DATAPOINTS);
        
        % Plot schematic picture of the Snowhite robot
        alfa = 660;
        beta = 0;
        gamma = -90*pi/180;
        figure(fig_env); plot_line_segments(REF, LINES, 1);
        plot_threewheeled_Laser([X(kk-1) Y(kk-1) AA(kk-1)]', 100, 612, 2, CONTROL(kk-1,5), 150, 50, 680, alfa, beta, gamma, angs, meas, 1);
        

        LINEMODEL = [REF(LINES(:,1),1:2) REF(LINES(:,2),1:2)];
 % Task 5a Kalman Filter with simulated position fixes (small) (Exercise 4)
        xv = 100;
        yv = 100;
        asv = 3*pi/180;
        C = [xv^2 0 0;
           0 yv^2 0;
           0 0 asv^2];
        x1 = randn;
        y1 = randn;
        a1 = randn;

        x2 = 0 + xv*x1;
        y2 = 0 + yv*y1;
        a2 = 0 + asv*a1;
        X_1 = CONTROL(kk-1,4) + x2;
        Y_1 = CONTROL(kk-1,5) + y2;
        A_1 = mod(CONTROL(kk-1,6) + a2, 2*pi);
        P1 = P(kk-1,1:9);      
    
        ICP = inv(C + reshape(P1,[3,3]));

        XYZ = [(X_1 -X(kk-1));
            (Y_1 -Y(kk-1));
            (AngDifference(A_1, AA(kk-1)))];
        
        XYZ__ = (reshape(P1,[3,3])*ICP)*XYZ;
        
        
        X(kk-1) = X(kk-1) + XYZ__(1,1);
        Y(kk-1) = Y(kk-1) + XYZ__(2,1);
        AA(kk-1) = mod(AA(kk-1) + XYZ__(3,1), 2*pi);
        
        c_XYZ = inv(inv(C) + inv(reshape(P1,[3,3])));
        
        P(kk-1,1:9) = reshape(c_XYZ,[1,9]);

 

        
        % Next time use the next scan
        scan_idx = mod(scan_idx, max(size(LD_HEAD))) + 1;
    end;
    
    % Mark the estimated (dead reckoning) position
    figure(fig_path);
    hold on; plot(X(kk-1), Y(kk-1), 'b.'); hold off;
    % Mark the true (from the LaserWay system) position
    hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'k.'); hold off;
    
    % Estimate the new position (based on the control inputs) and new
    % uncertainty
    V = CONTROL(kk-1,2);
    A = CONTROL(kk-1,3);
    T = 0.050;
    L = 680;

    X(kk) = X(kk-1) + V*cos(A)*T*cos(AA(kk-1) + (V*sin(A)*T/(2*L)));
    Y(kk) = Y(kk-1) + V*cos(A)*T*sin(AA(kk-1) + (V*sin(A)*T/(2*L)));
    AA(kk) = mod(AA(kk-1) + V*sin(A)*T/L, 2*pi);
    
    % ALSO UPDATE THE UNCERTAINTY OF THE POSITION
    % Task 2 Dead Reckoning - Add you code here


    % Task 2 Dead Reckoning - Add you code here
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];      
    Axya = [1 0 -T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 1  T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 0 1];
    Au = [(T*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T^2*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))  (- T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A) - (T^2*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L))  (V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
        (T*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T^2*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))   ((T^2*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L) - T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A))  (V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
                                                                                        (T*sin(A))/(L)                                                                                (T*V*cos(A))/(L)                                                                                 (V*sin(A))/(L)];
    
    Cu = [sv*abs(V) 0 0;
        0 sa^2 0;
        0 0 st^2];    
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Axya*Cxya_old*Axya' + Au*Cu*Au';

    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];    
end

ERROR = [X' Y' AA'] - CONTROL(:,4:6);
ERROR(:,3) = AngDifference(AA',CONTROL(:,6));
ERROR = abs(ERROR);
figure,
subplot(3,1,1);
plot(ERROR(:,1),'b'); hold;
plot(sqrt(P(:,1)),'r'); % one sigma
plot(ScanPosIndx,sqrt(P(ScanPosIndx,1)),'k.'); % one sigma
title('Error X [mm] and uncertainty [std] (red)');

subplot(3,1,2);
plot(ERROR(:,2),'b'); hold;
plot(sqrt(P(:,5)),'r'); % one sigma
title('Error Y [mm] and uncertainty [std] (red)');

subplot(3,1,3);
plot(ERROR(:,3)*180/pi,'b'); hold;
plot(sqrt(P(:,9))*180/pi,'r'); % one sigma
title('Error A [degree] and uncertainty [std] (red)');
%%
%%
% Task 6
%%
%
% Main script that reads controller data and laser data

clear all;
close all;

% Co-ordinates of the ref. nodes
REF = [1920 9470;   % 01
       10012 8179;  % 02
       9770 7590;   % 03
       11405 7228;  % 04
       11275 6451;  % 05
       11628 6384.5;% 06
       11438 4948;  % 07
       8140 8274;   % 08
       8392 8486;   % 09
       3280 2750;   % 10
       7250 2085;   % 11
       9990 1620;   % 12
       7485 3225;   % 13
       9505 3893;   % 14
       9602 4278;   % 15
       10412 4150;  % 16
       4090 7920;   % 17
       8010 5290;   % 18
       8255 6099;   % 19
       7733 6151;   % 20
       7490 6136;   % 21
       7061 5420;   % 22
       7634 5342];  % 23

% LINE SEGMENT MODEL (Indeces in the REF-vector)
LINES = [1 8;       % L01
         9 2;       % L02
         2 3;       % L03
         3 4;       % L04
         4 5;       % L05
         5 6;       % L06
         6 7;       % L07
         17 10;     % L08
         10 12;     % L09
         11 13;     % L10
         12 16;     % L11
         16 15;     % L12
         15 14;     % L13
         19 21;     % L14
         22 18;     % L15
         20 23];    % L16
         
% Control inputs (velocity and steering angle)
CONTROL = load('control_joy.txt');

% Laser data
LD_HEAD   = load('laser_header.txt');
LD_ANGLES = load('laser_angles.txt');
LD_MEASUR = load('laser_measurements.txt');
LD_FLAGS  = load('laser_flags.txt');

[no_inputs co] = size(CONTROL);

% Robots initial position
X(1) = CONTROL(1,4);
Y(1) = CONTROL(1,5);
AA(1) = CONTROL(1,6);
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];
sv = .01;
sa = .2; %least decimal
st = .05/100; % 1 percentage
scan_idx = 1;
fig_path = figure;
fig_env = figure;
ScanPosIndx = [];

% Plot the line model
figure(fig_env); plot_line_segments(REF, LINES, 1);
     
for kk = 2:no_inputs,
    % Check if we should get a position fix, i.e. if the time stamp of the
    % next laser scan is the same as the time stamp of the control input
    % values
    if LD_HEAD(scan_idx,1) == CONTROL(kk-1,1),
        % Mark the position where the position fix is done - and the size
        % of the position fix to be found
        figure(fig_path);
        hold on; plot(X(kk-1), Y(kk-1), 'ro'); hold off;
        hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'ro'); hold off;
        hold on; plot([X(kk-1) CONTROL(kk-1,4)], [Y(kk-1) CONTROL(kk-1,5)], 'r'); hold off;
        
        % > Log pos (index) there scan is taken
        ScanPosIndx = [ScanPosIndx kk-1]; % Used for plotting
        
        % Get the position fix - Use data points that are ok, i.e. with
        % flags = 0
        DATAPOINTS = find(LD_FLAGS(scan_idx,:) == 0);
        angs = LD_ANGLES(scan_idx, DATAPOINTS);
        meas = LD_MEASUR(scan_idx, DATAPOINTS);
        
        % Plot schematic picture of the Snowhite robot
        alfa = 660;
        beta = 0;
        gamma = -90*pi/180;
        figure(fig_env); plot_line_segments(REF, LINES, 1);
        plot_threewheeled_Laser([X(kk-1) Y(kk-1) AA(kk-1)]', 100, 612, 2, CONTROL(kk-1,5), 150, 50, 680, alfa, beta, gamma, angs, meas, 1);
        
        % Task 3 - write your code in the function "Cox_LineFit"
        % [dx,dy,da,C] = Cox_LineFit(angs, meas, [X(kk-1) Y(kk-1) A(kk-1)]', LINEMODEL,[alfa beta gamma]);
        % Returns => Position fix + Unceratinty of the position fix
        LINEMODEL = [REF(LINES(:,1),1:2) REF(LINES(:,2),1:2)];
        [dx,dy,da,C] = Cox_LineFit1(angs, meas, [X(kk-1) Y(kk-1) AA(kk-1)]', LINEMODEL,[alfa beta gamma]);
       
        % Task 6 Kalmanfilterr with Cox position update (Exercise 4)
        X_1 = X(kk-1) + dx;
        Y_1 = Y(kk-1) + dy;
        A_1 = mod(AA(kk-1) + da, 2*pi);
        P1 = P(kk-1,1:9);      


        ICP = inv(C + reshape(P1,[3,3]));

         XYZ = [(X_1 -X(kk-1));
            (Y_1 -Y(kk-1));
            (AngDifference(A_1, AA(kk-1)))];
        
        XYZ__ = (reshape(P1,[3,3])*ICP)*XYZ;
        
        X(kk-1) = X(kk-1) + XYZ__(1,1);
        Y(kk-1) = Y(kk-1) + XYZ__(2,1);
        AA(kk-1) = mod(AA(kk-1) + XYZ__(3,1), 2*pi);
        
        c_XYZ = inv(inv(C) + inv(reshape(P1,[3,3])));
        P(kk-1,1:9) = reshape(c_XYZ,[1,9]);

 

        
        % Next time use the next scan
        scan_idx = mod(scan_idx, max(size(LD_HEAD))) + 1;
    end;
    
    % Mark the estimated (dead reckoning) position
    figure(fig_path);
    hold on; plot(X(kk-1), Y(kk-1), 'b.'); hold off;
    % Mark the true (from the LaserWay system) position
    hold on; plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'k.'); hold off;
    
    % Estimate the new position (based on the control inputs) and new
    % uncertainty
    V = CONTROL(kk-1,2);
    A = CONTROL(kk-1,3);
    T = 0.050;
    L = 680;

    X(kk) = X(kk-1) + V*cos(A)*T*cos(AA(kk-1) + (V*sin(A)*T/(2*L)));
    Y(kk) = Y(kk-1) + V*cos(A)*T*sin(AA(kk-1) + (V*sin(A)*T/(2*L)));
    AA(kk) = mod(AA(kk-1) + V*sin(A)*T/L, 2*pi);
    
    % ALSO UPDATE THE UNCERTAINTY OF THE POSITION
    % Task 2 Dead Reckoning - Add you code here


    % Task 2 Dead Reckoning - Add you code here
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];      
    Axya = [1 0 -T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 1  T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A);
        0 0 1];
    Au = [(T*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T^2*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))  (- T*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A) - (T^2*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L))  (V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) - (T*V^2*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
        (T*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T^2*V*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L))   ((T^2*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)^2)/(2*L) - T*V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*sin(A))  (V*sin(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A) + (T*V^2*cos(AA(kk-1) + (T*V*sin(A))/(2*L))*cos(A)*sin(A))/(2*L));
                                                                                        (T*sin(A))/(L)                                                                                (T*V*cos(A))/(L)                                                                                 (V*sin(A))/(L)];
    
    Cu = [sv*abs(V) 0 0;
        0 sa^2 0;
        0 0 st^2];    
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Axya*Cxya_old*Axya' + Au*Cu*Au';

    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];    
end

ERROR = [X' Y' AA'] - CONTROL(:,4:6);
ERROR(:,3) = AngDifference(AA',CONTROL(:,6));
ERROR = abs(ERROR);
figure,
subplot(3,1,1);
plot(ERROR(:,1),'b'); hold;
plot(sqrt(P(:,1)),'r'); % one sigma
plot(ScanPosIndx,sqrt(P(ScanPosIndx,1)),'k.'); % one sigma
title('Error X [mm] and uncertainty [std] (red)');

subplot(3,1,2);
plot(ERROR(:,2),'b'); hold;
plot(sqrt(P(:,5)),'r'); % one sigma
title('Error Y [mm] and uncertainty [std] (red)');

subplot(3,1,3);
plot(ERROR(:,3)*180/pi,'b'); hold;
plot(sqrt(P(:,9))*180/pi,'r'); % one sigma
title('Error A [degree] and uncertainty [std] (red)');
%%
%%


