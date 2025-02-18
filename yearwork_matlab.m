load("structure_yw_mkr.mat");
%%
ndof = 68;
nvar = 72;
KFF = K(1:ndof,1:ndof);
KFC = K(1:ndof,ndof+1:nvar);
KCF = K(ndof+1:nvar,1:ndof);

MFF = M(1:ndof,1:ndof);
MFC = M(1:ndof,ndof+1:nvar);
MCF = M(ndof+1:nvar,1:ndof);

alpha = 0.2;
beta = 4*10^(-4);
CFF = R(1:ndof,1:ndof);
CFC = R(1:ndof,ndof+1:nvar);
CCF = R(ndof+1:nvar,1:ndof);


[eigenvect, eigenval] = eig(MFF\KFF);
modal_shapes = eigenvect;                       
modal_shapes_10hz = modal_shapes(:,55:60);       
modal_shapes_4 = modal_shapes(:,55:58);

modal_mat = modal_shapes;
modal_MFF = modal_mat' * MFF * modal_mat;
modal_KFF = modal_mat' * KFF * modal_mat;
modal_CFF = modal_mat' * CFF * modal_mat;

modal_mat = modal_shapes_4;
modal_MFF_4 = modal_mat' * MFF * modal_mat;       % modal mass matrix
modal_KFF_4 = modal_mat' * KFF * modal_mat;       % modal stiffness matrix
modal_CFF_4 = modal_mat' * CFF * modal_mat;       % modal damping matrix


%% 
% compute the natural frequencies of the damped structure up to 10 hz and
% the related non dimensional damping ratios
% Undamped case
omega_undamped = sqrt(diag(eigenval));
natural_freq_undamped = omega_undamped/pi/2;

% Damped case
non_dim_damp_ratio = diag(modal_CFF)/2./(diag(modal_MFF).*diag(modal_KFF)).^0.5;
omega_damped = natural_freq_undamped.*(1-non_dim_damp_ratio.^2).^0.5;   %omega_d = omega*sqrt(1-h^2)
natural_freq_damped_10hz = omega_damped/2/pi;

% Damped case up to 10hz
non_dim_damp_ratio_10hz = non_dim_damp_ratio(55:60);
natural_freq_undamped_10hz = natural_freq_undamped(55:60);
natural_freq_damped_10hz = natural_freq_damped_10hz(55:60);


%% 
% For the unloaded crane (i.e. without mass MA), develop a model in modal coordi-
% nates limited to the first four modes and plot the Bode diagrams (in linear scales)
% of the following frequency response functions (FRF) in the frequency range 0 ÷ 10
% Hz with step 0.01 Hz:
% Input: vertical force at point A; output: vertical displacement of point A;
% Input: horizontal force at point A; output: horizontal displacement of point
% B;
% Input: vertical force at point A; output: vertical acceleration of point A;

vett_f = 0:0.01:10;  % Frequency range

% Frequency Response Function (FRF) calculation for reduced system
mod1 = zeros(length(vett_f), 1);  
phase1 = zeros(length(vett_f), 1); 
mod2 = zeros(length(vett_f), 1);
phase2 = zeros(length(vett_f), 1);
mod3 = zeros(length(vett_f), 1);
phase3 = zeros(length(vett_f), 1);
i = sqrt(-1);

% Definition of forces
% a
F_a = zeros(ndof,1);              
F_a(20,1) = 1;                    
F_mod_a = modal_mat' * F_a; 

% b
F_b = zeros(ndof,1);              
F_b(19,1) = 1;                    
F_mod_b = modal_mat' * F_b; 

% c
F_c = zeros(ndof,1);              
F_c(20,1) = 1;                    
F_mod_c = modal_mat' * F_c; 

%% FRF1
for k = 1:length(vett_f)
    ome = vett_f(k) * 2 * pi;   
    A = -ome^2 * modal_MFF_4 + i * ome * modal_CFF_4 + modal_KFF_4; 
    q = A\ F_mod_a;  
    x = modal_mat*q;  
    
    phase1(k) = angle(x(20));   
    mod1(k) = abs(x(20));       
end

% Plot Magnitude Response
figure
subplot(2, 1, 1);
plot(vett_f, mod1, 'r', 'LineWidth', 0.7);
grid on;
xlabel('Frequency [Hz]');ylabel('Magnitude [N/m]');title('F_{Y_A}/Y_A')

% Plot Phase Response
subplot(2, 1, 2); hold on;
plot(vett_f, phase1 * 180 / pi, 'r', 'LineWidth', 0.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');


%% FRF2
for k=1:length(vett_f)
    ome=vett_f(k)*2*pi;
    A = -ome^2 * modal_MFF_4 + i * ome * modal_CFF_4 + modal_KFF_4; 
    q = A\F_mod_b;
    x = modal_mat*q;

    phase2(k)=angle(x(4));
    mod2(k)=abs(x(4));
end
figure
% Plot Magnitude Response
subplot(2, 1, 1);
plot(vett_f, mod2, 'r', 'LineWidth', 0.7);
grid on;
xlabel('Frequency [Hz]');ylabel('Magnitude [N/m]');title('F_{X_A}/X_B')

% Plot Phase Response
subplot(2, 1, 2);
plot(vett_f, phase2 * 180 / pi, 'r', 'LineWidth', 0.7); 
grid on;
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');

%% FRF3
for k=1:length(vett_f)
    ome=vett_f(k)*2*pi;
    A = -ome^2 * modal_MFF_4 + i * ome * modal_CFF_4 + modal_KFF_4; 
    q = A\F_mod_c;
    x = modal_mat*q;

    out = -ome^2*x(20);
    phase3(k)=angle(out);
    mod3(k)=abs(out);
end
figure
% Plot Magnitude Response
subplot(2, 1, 1);
plot(vett_f, mod3, 'r', 'LineWidth', 0.7); 
grid on;
xlabel('Frequency [Hz]');ylabel('Magnitude [N/(m/s^2]');title('F_{Y_A}/a_{Y_A}')

% Plot Phase Response
subplot(2, 1, 2); hold on;
plot(vett_f, phase3 * 180 / pi, 'r', 'LineWidth', 0.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');


%% 
% For the unloaded crane (i.e. without mass MA), plot the Bode diagrams (in linear scales) of the
% following frequency response functions (FRF) in the frequency range 0 ÷ 10 Hz with step 0.01
% Hz.
% a. Input: vertical force at point A; output: vertical component of the constraint force in
% the hinge O1;
% b. Input: horizontal force at point A; output: vertical component of the constraint force
% in the hinge O2.

i=sqrt(-1);
vett_f=0:0.01:10;
mod1 = zeros(length(vett_f), 1);
phase1 = zeros(length(vett_f), 1);
mod2 = zeros(length(vett_f), 1);
phase2 = zeros(length(vett_f), 1);

% Definition of forces
F_a = zeros(ndof,1);
F_a(20,1) = 1;
F_b = zeros(ndof,1);
F_b(19,1) = 1;

%  FRF 1 
for k=1:length(vett_f)
    ome=vett_f(k)*2*pi;
    A=-ome^2*MFF+i*ome*CFF+KFF;
    x = A\F_a;
    F_co = (-ome^2*MCF+i*ome*CCF+KCF)*x;
    phase1(k)=angle(F_co(2));
    mod1(k)=abs(F_co(2));
end
figure
subplot 211; plot(vett_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Magnitude [N/N]');title('F_{Y_A}/Fc_{yO_1}')
subplot 212; plot(vett_f,phase1*180/pi);grid;xlabel('Frequency [Hz]');ylabel('Phase [deg]')


for k=1:length(vett_f)
    ome=vett_f(k)*2*pi;
    A=-ome^2*MFF+i*ome*CFF+KFF;
    x = A\F_b;
    F_co = (-ome^2*MCF+i*ome*CCF+KCF)*x;
    phase1(k)=angle(F_co(4));
    mod1(k)=abs(F_co(4));
end
figure
subplot 211; plot(vett_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Magnitude [N/N]');title('F_{X_A}/Fc_{yO_2}')
subplot 212; plot(vett_f,phase1*180/pi);grid;xlabel('Frequency [Hz]');ylabel('Phase [deg]')





%% 
% For the unloaded crane (i.e. without mass MA), plot the Bode diagram (in linear scales) of the
% following frequency response function (FRF) in the frequency range 0 ÷ 10 Hz with step 0.01
% Hz.
% a. Input: horizontal force at point A; output: bending moment at section C-C, belonging
% to the vertical green beam and located just above the connection with the red
% transverse beam

mod1 = zeros(length(vett_f), 1);
phase1 = zeros(length(vett_f), 1);
i=sqrt(-1);
vett_f=0:0.01:10;

F = zeros(ndof,1);
F(19,1) = 1;

E = 2.06*10^11;
I = 6.712*10^-4;
EI = E*I;
L = 6;
for k = 1:length(vett_f)
    ome = vett_f(k) * 2 * pi; 
    A = -ome^2*MFF+i*ome*CFF+KFF;
    x = A \ F;

    %  Find displacement and rotation of boundaries of C-C section (nodes
    %  8 and 23)
    % Nodal coordinates
    xgi = x(22); ygi = x(23); thetagi = x(24);
    xgj = x(63); ygj = x(64); thetagj = x(65);

    % Transformation to local reference
    xiL = ygi; yiL = -xgi; thetaiL = thetagi;
    xjL = ygj; yjL = -xgj; thetajL = thetagj;

    % solving the system of boundaries condition using cubic
    % interpolating function
    % w(x) = ax^3+bx^2+cx^+d -> vertical displacement
    % w'(x) = 3ax^2+2bx+c -> angular displacement
    % w''(x) = 6ax+2b -> bending moment = EI*w''(x)

    % Boundary conditions
    % w(0) = yiL | w(L) = yjL | w'(0) = thetaiL | w'(L) = thetajL

    boundaries_matrix = [0      0       0       1;
                         0      0       1       0;
                         L^3    L^2     L       1;
                         3*L^2  2*L     1       0];

    local_coordinates = [yiL;  thetaiL;  yjL;  thetajL];

    coefficients = boundaries_matrix\local_coordinates;
    
    M_CC = EI*2*coefficients(2);

    % Store magnitude and phase
    mod1(k) = abs(M_CC);
    phase1(k) = angle(M_CC);
end

figure
subplot(2, 1, 1);
plot(vett_f, mod1); grid;
xlabel('Frequency [Hz]'); ylabel('Magnitude [N/Nm]');
title('Bending moment section C-C');

subplot(2, 1, 2);
plot(vett_f, phase1*180/pi); grid;
xlabel('Frequency [Hz]'); ylabel('Phase [deg]');


