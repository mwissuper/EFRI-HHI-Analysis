% Manual check on recoverForces.m code
% For Fx (vicon CS), check for rotation of sensor about y and z axes and
% take projections of component forces on x (vicon) axis
% Use front left and front right FH markers
% First run mainWorkPowerAnalysisMW.m to just after first call of
% recoverForces.m

% Scale Factors (Convert V to N or Nm)
F_Scale = [12.5, 12.5, 50];     %Coordinates are scaling factors in [X, Y, Z]
T_Scale = [0.3, 0.3, 0.3];
% Sensor Biases (Sensor offsets in V)
F_Bias = [-0.2567, -0.1931, -0.1590]; %Force bias in [X, Y, Z] (V) coordinates
T_Bias = [-0.6627, 0.6189, 0.7339]; %Torque bias
% Convert sensor output (V) to forces (N) and torques (Nm) assuming a
% linear relationship: (Force/Torque) = (Scale)[(Sensor Reading) -
% (Bias)]. This conversion is performed in the Force/Torque Sensor's
% Reference frame
F = F_Scale.*(Force - F_Bias);
T = T_Scale.*(Torque - T_Bias);

% Get rotation from markers on handle
a = FH.frontright' - FH.frontleft';

rotz = atan2(a(:,2),a(:,1));
roty = atan2(a(:,3),a(:,1));
% when angle of rotation is small, the x component of the transducer's
% force is very similar to the recovered force's x component.

subplot(2,1,1)
plot(rotz),hold on, plot(roty); % Plot to see when rotation angles are large or small
legend('z','y'); ylabel('Rotation angle (rad)');
B = [1 0 0]; % Lab CS x dir
for i = 1:length(FH.frontright') % for each time point of marker data
    A = F(i,:); % Look at 3D force vector at this time point
    Fx(i) = dot(A,B); % Just pull out magnitude, already know what dir vector is in
end

subplot(2,1,2)
plot(Fx),hold on;

% Use recoverForces code to compare calculation of Fx
ForceRF = recoverForces(Force',Torque',FH); 
ForceRF = ForceRF';

plot(ForceRF(:,1))
legend('Fx','Fx rF');