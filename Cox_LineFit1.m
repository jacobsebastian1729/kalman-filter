%
% Fits (by translation and rotation) data points to a set of 
% line segments, i.e. applying the Cox's algorithm
% function [ddx,ddy,dda,C] = Cox_LineFit(ANG, DIS, POSE, LINEMODEL, SensorPose)
% 
function [ddx,ddy,dda,C] = Cox_LineFit1(ANG, DIS, POSE, LINEMODEL, SensorPose)
    % Init variables
    ddx = 0; ddy = 0; dda = 0;
    Rx = POSE(1,1); Ry = POSE(2,1); Ra = POSE(3,1); 
    sALFA = SensorPose(1); sBETA = SensorPose(2); sGAMMA = SensorPose(3);
    max_iterations = 100; % <------YOU NEED TO CHANGE THIS NUMBER
    no_update = 0;
    
    % Step 0 - Normal vectors (length = 1) to the line segments
    % -> Add your code here
    
    a1 = LINEMODEL(:,3) -LINEMODEL(:,1);
    b1 = LINEMODEL(:,4) -LINEMODEL(:,2);
    c1 = sqrt(a1.^2 + b1.^2);
    uv(:,1) = -b1./c1;
    uv(:,2) = a1./c1;
    tr = (LINEMODEL(:,1:2));
    ri = dot(transpose(tr),transpose(uv));
    % The Loop - REPEAT UNTIL THE PROCESS CONVERGE
    %b(1) = 0;b(2) = 0;b(3) = 0;
    for iteration = 1:100,
        % Step 1 Translate and rotate data points
            % 1.1) Relative measurements => Sensor co-ordinates
            %-> Add your code here
            x1 = DIS.*cos(ANG);
            y1 = DIS.*sin(ANG);
            % 1.2) Sensor co-ordinates => Robot co-ordinates
            %-> Add your code here
            R1 = [cos(sGAMMA) -sin(sGAMMA) sALFA; sin(sGAMMA) cos(sGAMMA) sBETA;0 0 1];
            z1 = ones(1, length(x1));
            xyz1(1,:) = x1;
            xyz1(2,:) = y1;
            xyz1(3,:) = z1;
            Xs = R1*xyz1;
            % 1.3) Robot co-ordinates => World co-ordinates
            %-> Add your code here
            R2 = [cos(Ra + dda) -sin(Ra + dda) Rx+ddx; sin(Ra + dda) cos(Ra + dda) Ry+ddy;0 0 1];
            Xs_ = Xs;
            Xw_ = R2*Xs_;
            Xw = Xw_;
            Xw1 = Xw(1:2,:);
            yi = transpose(ri) - uv*Xw1;
            %disp(Xw);
            %plot(Xw1(1,:), Xw1(2,:));
        % Step 2 Find targets for data points
        %-> Add your code here
            yi2 = abs(yi);
            minyi(1,:) = min(yi2);
            median_ = median(minyi);
            %disp(median_);
            for i = 1:length(minyi)
                index_r(1,i) = find(yi2(:,i) == min(yi2(:,i)));
                yi4(1,i) = yi(index_r(1,i),i);
            end    
    
            index_ = find(minyi <= median_);
            minY = minyi(index_);
            minYi = yi4(index_);
            ang1 = ANG(index_);
            mes1 = DIS(index_);
            for i = 1:length(yi2)
                yi3(1,i) = find(yi2(:,i) == min(yi2(:,i)));
            end
        % Step 3 Set up linear equation system, find b = (dx,dy,da)' from the LS
        %-> Add your code here
            xi1_ = uv(yi3,1);
            xi1 = xi1_(index_);
            xi2_ = uv(yi3,2);
            xi2 = xi2_(index_);
            uv_xi(:,1) = xi1;
            uv_xi(:,2) = xi2;
            cg(1,:) = sum(Xw(1,index_),2)/length(Xw(1,index_));
            cg(2,:) = sum(Xw(2,index_),2)/length(Xw(1,index_));
            xw_cg(1,:) = Xw(1,index_);
            xw_cg(2,:) = Xw(2,index_);
            for i = 1:length(Xw(1,index_))
                xi3(i,1) = uv_xi(i,:)*[0 -1;1 0]*(xw_cg(:,i)-cg);
            end

            A1 = [xi1 xi2 xi3];
            B1 = inv(transpose(A1)*A1);
            b = B1*transpose(A1)*transpose(minYi);
        %b = POSE/max_iterations; % <--- You shall change this! This is only
        % for demonstation, i.e. return the same pose as sent in.
            s_0 = minYi-transpose(A1*b);
            s_2 = s_0*(transpose(s_0)) / (max(size(A1)) - 4);

            C = s_2*inv(transpose(A1)*A1); % Covarince matrix
        
        % Step 4 Add latest contribution to the overall congruence 
            ddx = ddx + b(1);
            ddy = ddy + b(2);
            dda = mod(dda + b(3), 2*pi);
        % Step 5  Check if the process has converged
        %-> Add your code here
            if (sqrt(b(1)^2 + b(2)^2) < 5 )&&(abs(b(3) <0.1*pi/180))
                %disp(iteration)
                break
            end

    end;
   
    %disp(b);
    %disp(iteration)
    