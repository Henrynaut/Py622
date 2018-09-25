%Written by Neil McHenry
%September 23rd, 2018

%AERO 622 Hw %3 Page 49 Problem 2
    %Solve for the motion of the bead on a spinning wire
    %Using a numerical integration technique, then comment

clc; clear variables; close all; tic;

global OMEGA m g r;

%Physical parameters
    m = 0.1;                  %kg
    r = 0.4;                  %meters
    g = 10;                   %m/s/s
    thetaDot = 0;             %rad/s
    theta = -pi/2 + 0.1;      %radians
         
%ODE45 Parameters
    options=odeset('RelTol', 1e-5);   
    t0 = 0;
    initial = [theta, thetaDot];
  

%% First set of Conditions
    OMEGA = 4;                  %rad/s
    t1 = 12;                    %seconds
    time=[t0,t1];               %Array of start and end times
        
    %Solve with ODE45
    [A1,B1]=ode45(@beadWire,time,initial,options);    

%% Second set of Conditions
    OMEGA = 6;              %rad/s    
    t1 = 12;               %seconds
    time=[t0,t1];      

    %Solve with ODE45
    [A2,B2]=ode45(@beadWire,time,initial,options);     
    
    %% Plots

figure (1)
subplot(1,2,1)
hold on
plot(A1,rad2deg(B1(:,1)))
plot(A2,rad2deg(B2(:,1)), '*-')
xlabel('Time, s')
ylabel('Degrees')
title('Theta of Bead on Spinning Wire')
legend('OMEGA = 4 rad/s.', 'OMEGA = 6 rad/s.')
hold off

subplot(1,2,2)
hold on
plot(A1,rad2deg(B1(:,2)))
plot(A2,rad2deg(B2(:,2)), '*-')
xlabel('Time, s')
ylabel('Degrees per Second')
title('Theta Dot of Bead on Spinning Wire')
legend('OMEGA = 4 rad/s.', 'OMEGA = 6 rad/s.')
hold off

    
