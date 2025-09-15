%Function to calculate structural mass of electrical machines
%Example:
%For the 10 MW,10 rpm machine
%[rotor_mass, stator_mass]=structure(9.5e6,4.6,5.8,1.4,1.8)
%rotor mass = 42.8 t, stator = 82.6 t
%For the 36.5 MW, 120 rpm machine
%[rotor_mass, stator_mass]=structure(2.9e6,3.5,4,1.3,1.8)
%rotor mass = 29.8 t, stator = 48 t

%   INPUTS
% T: Mechanical Torque (Nm)
% D: Rotor Structure Diameter (m)
% Ds: Stator Inner Diameter (m)
% L: Axial Length of Machine (m)
% B_g: Maximum value of air gap flux density(T) (for normal stress calculation)

%    OUTPUTS
% mass_rotor_structure: Total structural mass of rotor (!Mass of shaft should be added)
% mass_stator_structure: Total structural mass of stator (!Mass of bearing may be added)

function [mass_rotor_structure, mass_stator_structure] = PM_structure(T, D, Ds, L, B_g)

% Inputs from Other modules
% T=2.3e6;     %Torque in Nm
% 
% D = 6;       % Rotor Inner Diameter (m) without magnet height
% Ds = 6.04;   % Stator Inner Diameter (m) -- (Dr + 2*magnet_height + 2*airgap) 
% L = 1.7;     % Axial Length (m)
% 
% B_g = 1.1;   % Maximum Flux Density in the Airgap (T)

%% Main Variables %%

%should be linked with number of poles and diameter
n = 5;          % Number of Structural Arms for the Rotor
n_s = 5;        % Number of Structural Arms for the Stator

%% Constants %%

E = 2e11;          % Young's Modulus (Pa)
rho_steel = 7800;  % Density of Steel (kg/m^3)
mu_0 = pi.*4e-7;      % Permeability of Free Space (H/m)

% Stresses %

q = (B_g.^2)./2/mu_0;   % Normal Component of Maxwell Stress (Pa)

sigma = 2*T./pi./L./Ds.^2;       % Sheer Stress (Pa)


%May change according to machine topology enter as an input?


R_o = D./10;         % Shaft Radius

hry = D./300;        % Rotor Yoke Height
hsy = Ds./300;       % Stator Yoke Height
b = D.*R_o/n;        % Hollow Rotor Arm's Width
d = D./200;          % Hollow Rotor Arm's Length 
t_w = D./200;        % Hollow Rotor Arm's Thickness
b_s = Ds.*R_o/n_s;   % Hollow Stator Arm's Width
d_s = Ds./200;       % Hollow Stator Arm's Length
t_ws = D./200;       % Hollow Stator Arm's Thickness


%% Deflections %%

% The deflections are needed as a criterion to ensure that the structures
% are strurdy enough against the magnetic attraction & torque

% Radial Deflection %

theta = pi./n;
theta_s = pi./n_s;
a = (b.*d)-((b-2.*t_w).*(d-2.*t_w));             % Hollow Rotor Arm Area (m^2)
a_s = (b_s.*d_s)-((b_s-2.*t_ws).*(d_s-2.*t_ws)); % Hollow Stator Arm Area (m^2)
I = L.*hry.^3/12;                                % Second Moment of Area of the Rotor�s Cylinder (m^4)
Is = L.*hsy.^3/12;                               % Second Moment of Area of the Stator�s Cylinder (m^4)
A = L.*hry;                                      % Cross Area of the Rotor�s Cylinder (m^2)
As = L.*hsy;                                     % Cross Area of the Stator�s Cylinder (m^2)
k = sqrt(I./A);                                  % Radius of Rotor Gyration (m)
mnt = (k./(D./2)).^2;
k_s = sqrt(Is./As);                              % Radius of Stator Gyration (m)
mnt_s = (k_s./(Ds./2)).^2;
u_r = (q.*(D./2).^2/E./hry).*(1+(((D./2).^3*((0.25.*(sin(theta)-(theta.*cos(theta)))./(sin(theta)).^2)...
    -(0.5./sin(theta))+(0.5./theta)))./I./((((theta./(sin(theta).^2))+(1./tan(theta))).*((0.25.*(D./2)./A)...
    +(0.25.*(D./2).^3/I)))-((D./2).^3/(2.*I.*theta.*(mnt+1)))+((D./2-hry-R_o)./a))));                                 % Radial Deflection of Rotor due to Maxwell Stress (m)
u_s = (q.*(Ds./2).^2/E./hsy).*(1+(((Ds./2).^3*((0.25.*(sin(theta_s)-(theta_s.*cos(theta_s)))./(sin(theta_s)).^2)...
    -(0.5./sin(theta_s))+(0.5./theta_s)))./Is./((((theta_s./(sin(theta_s).^2))+(1./tan(theta_s))).*((0.25.*(Ds./2)./As)...
    +(0.25.*(Ds./2).^3/Is)))-((Ds./2).^3/(2.*Is.*theta_s.*(mnt_s+1)))+((Ds./2-hsy-R_o)./a_s.*2))));                    % Radial Deflection of Stator due to Maxwell Stress (m)

u_r_all = D./2./1e4;          % Allowable Radial Deflection of the Rotor (m)
u_s_all = Ds./2./1e4;         % Allowable Radial Deflection of the Stator (m)

% Tangential Deflection %

I_arm_tor = ((d.*b.^3)-((d-2.*t_w).*(b-2.*t_w).^3))./12;                 % Second Moment of the Area of the Rotor's Arm
I_arm_tor_s = ((d_s.*b_s.^3)-((d_s-2.*t_ws).*(b_s-2.*t_ws).^3))./12;   % Second Moment of the Area of the Stator's Arm
l_i = (D./2) - R_o;                                              % Length of Rotor Arm (m)
l_i_s = (Ds./2) - R_o;                                           % Length of Stator Arm (m)
z = 2*pi*(D./2)*L.*sigma.*(l_i.^3)./n/3/E/I_arm_tor;              % Circumferential Deflection of the Rotor Arm (m)
z_s = (pi.*(Ds./2).*L./n_s).*sigma.*(l_i_s.^3)./3./E./I_arm_tor_s;        % Circumferential Deflection of the Stator Arm (m)

z_all = 0.005.*2.*pi.*(D./2)./360;       % Allowable Circumferential Deflection of the Rotor Arm (m)
z_s_all = 0.005.*2.*pi.*(Ds./2)./360;    % Allowable Circumferential Deflection of the Stator Arm (m)

while z_all < z
        d = d+0.002;
        I_arm_tor = ((d.*b.^3)-((d-2.*t_w).*(b-2.*t_w).^3))./12;
        l_i = (D./2) - R_o;
        z = 2.*pi.*(D./2).*L.*sigma.*(l_i.^3)./n./3./E./I_arm_tor;
end
 
while u_r_all < u_r
          hry = hry+0.002;
          a = (b.*d)-((b-2.*t_w).*(d-2.*t_w));
          I = L.*hry.^3/12;
          A = L.*hry;
          k = sqrt(I./A);
          mnt = (k./(D./2)).^2;
          u_r = (q.*(D./2).^2/E./hry).*(1+(((D./2).^3*((0.25.*(sin(theta)-(theta.*cos(theta)))./(sin(theta)).^2)...
                -(0.5./sin(theta))+(0.5./theta)))./I./((((theta./(sin(theta).^2))+(1./tan(theta))).*((0.25.*(D./2)./A)...
                +(0.25.*(D./2).^3/I)))-((D./2).^3/(2.*I.*theta.*(mnt+1)))+((D./2-hry-R_o)./a))));

end
    
while z_s_all < z_s
         d_s = d_s+0.002;
         I_arm_tor_s = ((d_s.*b_s.^3)-((d_s-2.*t_ws).*(b_s-2.*t_ws).^3))/12;
         l_i_s = (Ds./2) - R_o;
         z_s = (pi.*(Ds./2).*L./n_s).*sigma.*(l_i_s.^3)./3./E./I_arm_tor_s;
end

while u_s_all < u_s
        hsy = hsy+0.002;
        a_s = (b_s.*d_s)-((b_s-2.*t_ws).*(d_s-2.*t_ws));
        Is = L.*hsy.^3/12;
        As = L.*hsy;
        k_s = sqrt(Is./As);
        mnt_s = (k_s./(Ds./2)).^2;
        u_s = (q.*(Ds./2).^2/E./hsy).*(1+(((Ds./2).^3*((0.25.*(sin(theta_s)-(theta_s.*cos(theta_s)))./(sin(theta_s)).^2)...
              -(0.5./sin(theta_s))+(0.5./theta_s)))./Is./((((theta_s./(sin(theta_s).^2))+(1./tan(theta_s))).*((0.25.*(Ds./2)./As)...
              +(0.25.*(Ds./2).^3/Is)))-((Ds./2).^3/(2.*Is.*theta_s.*(mnt_s+1)))+((Ds./2-hsy-R_o)./a_s.*2))));
end
    

%% Outputs %%

%mass_rotor_structure = ((pi.*hry.*L.*D.*rho_steel)+(n.*(D./2-hry-R_o).*a.*rho_steel));        %Rotor mass excluding PM and steel laminations
mass_rotor_structure = (n.*(D./2-hry-R_o).*a.*rho_steel);        %Rotor mass excluding PM and steel laminations
mass_stator_structure = (n_s.*(Ds./2-hsy-R_o).*2.*a_s.*rho_steel); %Stator mass excluding Copper and steel laminations

%total_mass = mass_rotor_structure + mass_stator_structure; 

%% Print Results %%
%fprintf('MASS:\nMass of rotor structure = %d kg\nMass of stator structure = %d kg\nTotal Mass = %d kg\n\n',round(mass_rotor_structure), round(mass_stator_structure), round(total_mass))

end



