%{
*
 NAME           : DiscreteKalmanFilter.m
 AUTHOR         : Paweekorn Buasakorn
 DATE           : May 07th 2023
 MODIFIED BY    : Paweekorn Buasakorn
 DESCRIPTION    : This function use to estimate Angular velocity and
                  current of Brushed DCmotor By the input is control input
                  and current position By Given State as [x1,x2,x3,x4] = [q,qd,m,i]
*
%}
global X_pk;
global P_pk;

[A,B,C,D,G] = MatrixGenerator();
X_pk = zeros(4,1);
P_pk = 0;

% 
% J = 6.38593361359858E-4;
% b = 0.0020388649255862683;
% kt = 0.11589677272727274;
% ke = 0.239;
% R = 1.119;
% L = 0.00466118;

% J = 1.7E-6;
% kt = 34.8E-3;
% ke = (3.65e-3)/60;
% R = 2.84;
% L = 380E-6;
% b = 4.37E-4;

%Parameter ZGB102FGG
J = 0.007371752058300;
kt =  0.395367525411770;
ke =  0.395367525411770;
R = 1.325581400000000;
L = 0.002386100000000;
b = 0.013503282026890;

G1 = 0.9;
G2 = 1- G1;
for i = 1:length(t)
    [w_hat(i),i_hat(i),p_hat(i)] = SteadyStateKalmanFilter(Vin(i),pos(i));
    i_hat(i) = G1*i_hat(i) + G2*(i_hat(i) - Current(i));
end
subplot(2,3,1);
plot(t,w_hat);
title('SteadyStateKalmanFilter [Motor model]')
xlabel('time (sec)')
ylabel('prediction velocity')

subplot(2,3,3);
plot(t,w);
xlabel('time (sec)')
ylabel('measurement velocity')

subplot(2,3,2);
plot(t,i_hat);
title('SteadyStateKalmanFilter [Motor model]')
xlabel('time (sec)')
ylabel('prediction current')

subplot(2,3,4);
plot(t,Current);
xlabel('time (sec)')
ylabel('measurement current')

subplot(2,3,5);
plot(t,p_hat);
subplot(2,3,6);
plot(t,pos);
function [A,B,C,D,G] = MatrixGenerator()
    % J = 922.333e-6;
    % b = 208.534e-6;
    % kt = 0.171;
    % ke = 0.239;
    % R = 1.119;
    % L = 5.1018e-3;
    

    %Parameter ZGB102FGG
    J = 0.007371752058300;
    kt =  0.395367525411770;
    ke =  0.395367525411770;
    R = 1.325581400000000;
    L = 0.002386100000000;
    b = 0.013503282026890;

    dt = 1.0/1000.0;

    %State Transition Matrix
    Ac = [0   1    0    0;
         0 -b/J -1/J kt/J;
         0   0    0    0;
         0 -ke/L  0  -R/L];
    %Input Matrix
    Bc = [0;
         0;
         0;
         1/L]; 
    %Process noise
    Gc = [0;
         1;
         0;
         0];
    %Output matrix
    Cc = [1 0 0 0];
    %Feed through matrix
    Dc = 0;
    
    % Create state space model
    sys = ss(Ac,Bc,Cc,Dc);
    % Convert to discrete, where dt is your discrete time-step (in seconds)
    d_sys = c2d(sys,dt);
    
    %Discrete Matrix
    A = d_sys.A;
    B = d_sys.B;
    C = d_sys.C;
    D = d_sys.D;
    G = Gc;
end


function [Es_Velocity ,Es_Current,Es_Position] = SteadyStateKalmanFilter(Vin,Position)
    global X_pk
    global P_pk
    Q = 1;
    R = 1;
    
    [A,B,C,D,G] = MatrixGenerator();
    %Prediction of the state Covariences error
    X_k = (A * X_pk) + (B * Vin);
    P_k = (A * P_pk * A.') + G*Q*G';
    
    %Computation of the Kalman Gain
    K = (P_k*C.')/(C*P_k*C.' + R);

    %Computation of estimative
    X_k = X_k + K*(Position - C*X_k);
    Es_Position = X_k(1);
    Es_Velocity = X_k(2);
    Es_Current = X_k(4);
    %Computation Covarrience error
    P_k = (eye(4)-K*C)*P_k; 

    X_pk = X_k;
    P_pk = P_k;
end
