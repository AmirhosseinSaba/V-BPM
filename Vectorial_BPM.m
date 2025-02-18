%% Vectorial Beam Propagation Method
%%%%%%%%%%%%%% Vectorial BPM %%%%%%%%%%%%%%

% Written By A. Saba, EPFL
% 20.05.2020
%
% This code presents vectorial beam propagation method for to calcualte
% Jones matrix and vectorial scattering from a inhomogenous anisotropic
% material. The anisotropic material is represented as a phantom with 9
% elements of its refractive index tensor and each element is a 3D array of Nz*Ny*Nx pixels. 
%
% Theoretical details in [1]
%
% Cite [1] in order to use this code
% [1] Saba, Amirhossein, et al. "Polarization-sensitive optical diffraction tomography." Optica 8.3 (2021): 402-408.


%% Parameters explaination

% dx: pixel size in um
% Nx, Ny, Nz: Number of pixels in X, Y, and Z direction
% wl: wavelength in um

%%

clear ; clc;close all;

Data_path           = pwd;
Data_path_s         = [Data_path, '\Phantom\'];
% Data_path_s         = [Data_path,'\'];

Sample_name         = 'Phantom_11_4_2020';

load([Data_path_s, Sample_name, '.mat']);

%%
[X,Y,Z] = meshgrid( -round((Nx+1)/2) + 1 : Nx - round((Nx+1)/2), -round((Ny+1)/2) + 1 : Ny - round((Ny+1)/2), -round((Nz+1)/2) + 1 : Nz - round((Nz+1)/2));                        % Center of phantom in the dimension of [y,x,z] in [um]

%% Defining RI and Physical Parameters %%
n0 = 1.333;
k0 = 2*pi./wl;
k  = k0*n0;
NA = 1;
effective_NA = min(NA,n0);
dz = dx;

%% Phantom Permittivity Tensor  (Not anymore be used in the code with the approximations)

epsilon11 = 2*n0*n11+n11.*n11+n12.*n21+n13.*n31;
epsilon12 = 2*n0*n12+n11.*n12+n12.*n22+n13.*n32;
epsilon13 = 2*n0*n13+n11.*n13+n12.*n23+n13.*n33;
epsilon21 = 2*n0*n21+n21.*n11+n22.*n21+n23.*n31;
epsilon22 = 2*n0*n22+n21.*n12+n22.*n22+n23.*n32;
epsilon23 = 2*n0*n23+n21.*n13+n22.*n23+n23.*n33;
epsilon31 = 2*n0*n31+n31.*n11+n32.*n21+n33.*n31;
epsilon32 = 2*n0*n32+n31.*n12+n32.*n22+n33.*n32;
epsilon33 = 2*n0*n33+n31.*n13+n32.*n23+n33.*n33;

%% Spatial frequencies

Kx          = 2*pi/Nx/dx*X(:,:,1);
Ky          = 2*pi/Ny/dx*Y(:,:,1);
Kx          = (Kx);
Ky          = (Ky);
K2          = Kx.^2 + Ky.^2;
K2          = ifftshift(K2);
Kz          = real(sqrt(k^2 - K2));
DFRa        =  zeros(Nx,Ny,Nz,2,2);

Numerical_Aperture     = (Kx.^2 + Ky.^2) < (k0*effective_NA)^2;     % Numerical Aperture of the System
Numerical_Aperture     = fftshift(Numerical_Aperture);
%%
DFR         = exp(-(1i) * K2 * dz ./ (k + Kz));                     % Diffraction Kernel
DFR_Bpr     = conj(DFR);                                            % Kernel of one pixel Backpropagation in z-axis
DFR_Mea     = exp(-(1i) * K2 * dz * floor((Nz-1)/2) ./ (k + Kz));	% Kernel of backpropagation to the center point in z-axis
DFR_Ori     = conj(exp(-(1i) * K2 * dz * round((Nz-1)/2) ./ (k + Kz)));
DFR_Ori2     = conj(exp(-(1i) * K2 * dz * round((Nz+1)/2) ./ (k + Kz)));

DFR         = (single(DFR));
DFR_Bpr     = (single(DFR_Bpr));
DFR_Mea     = (single(DFR_Mea));
DFR_Ori     = (single(DFR_Ori));

%%
theta = 25/180*pi;
eta = linspace(0,356,90);
eta = pi/180*eta;
Nview = length(eta);
kx_in = k*sin(theta).*cos(eta);
ky_in = k*sin(theta).*sin(eta);
kz_in = sqrt((k0*n0).^2 - kx_in.^2 - ky_in.^2);

%% Output Rotational Matrix
% Rout = Ny*Nx*3*3 Matrix, each one is indexed with Rij
% We generate each of these elements, Rij, which is a 2D array

Kz = real(sqrt(k^2-Kx.^2-Ky.^2));
thetaa = acos(Kz./(k));
ux = Kx/k;
uy = Ky/k;
uz = Kz/k;
u  = cat(3,ux,uy,uz);                                               % unitary vector along the propagation direction of the corresponding plane-wave
% calculation of the cross-product of zhat and u
v1 = -uy;                                                           
v2 = ux;
v3 = zeros(size(uz));
%
Rout11 = cos(thetaa)+ (v1.^2).*(1-cos(thetaa));
Rout12 = v1.*v2.*(1-cos(thetaa));
Rout13 = v2.*sin(thetaa);

Rout21 = v1.*v2.*(1-cos(thetaa));
Rout22 = cos(thetaa)+ (v2.^2).*(1-cos(thetaa));
Rout23 = -v1.*sin(thetaa);

Rout31 = -v2.*sin(thetaa);
Rout32 = v1.*sin(thetaa);
Rout33 = cos(thetaa);

Rout11 = (fftshift(Rout11));
Rout12 = (fftshift(Rout12));
Rout13 = (fftshift(Rout13));
Rout21 = (fftshift(Rout21));
Rout22 = (fftshift(Rout22));
Rout23 = (fftshift(Rout23));
Rout31 = (fftshift(Rout31));
Rout32 = (fftshift(Rout32));
Rout33 = (fftshift(Rout33));

Routinv11 = Rout11;
Routinv12 = Rout21;
Routinv13 = Rout31;
Routinv21 = Rout12;
Routinv22 = Rout22;
Routinv23 = Rout32;
Routinv31 = Rout13;
Routinv32 = Rout23;
Routinv33 = Rout33;

%% Input Rotational Matrix
% Rout = 3*3*TotalAng Matrix, each one is indexed with Rij
% We generate each of these elements, Rij, for each angle

uinx = kx_in/k;
uiny = ky_in/k;
uinz = kz_in/k;

theta_in = acos(uinz);

% calculation of the cross-product of zhat and u
vin1 = -uiny./(sqrt(uiny.^2+uinx.^2));                                                           
vin2 = uinx./(sqrt(uiny.^2+uinx.^2));
vin3 = zeros(size(uinz));
%
Rin11 = cos(theta_in)+ (vin1.^2).*(1-cos(theta_in));
Rin12 = vin1.*vin2.*(1-cos(theta_in));
Rin13 = vin2.*sin(theta_in);

Rin21 = vin1.*vin2.*(1-cos(theta_in));
Rin22 = cos(theta_in)+ (vin2.^2).*(1-cos(theta_in));
Rin23 = -vin1.*sin(theta_in);

Rin31 = -vin2.*sin(theta_in);
Rin32 = vin1.*sin(theta_in);
Rin33 = cos(theta_in);

Rin_inv11 = Rin11;
Rin_inv12 = Rin21;
Rin_inv13 = Rin31;
Rin_inv21 = Rin12;
Rin_inv22 = Rin22;
Rin_inv23 = Rin32;
Rin_inv31 = Rin13;
Rin_inv32 = Rin23;
Rin_inv33 = Rin33;
%%
cos_factor            = cos(theta);
tic
PHMa(:,:,:,1,1)       = 1i*k0*n11*dz/cos_factor;
PHMa(:,:,:,1,2)       = 1i*k0*n12*dz/cos_factor;
PHMa(:,:,:,1,3)       = 1i*k0*n13*dz/cos_factor;
PHMa(:,:,:,2,1)       = 1i*k0*n21*dz/cos_factor;
PHMa(:,:,:,2,2)       = 1i*k0*n22*dz/cos_factor;
PHMa(:,:,:,2,3)       = 1i*k0*n23*dz/cos_factor;
PHMa(:,:,:,3,1)       = 1i*k0*n31*dz/cos_factor;
PHMa(:,:,:,3,2)       = 1i*k0*n32*dz/cos_factor;
PHMa(:,:,:,3,3)       = 1i*k0*n33*dz/cos_factor;

PHM                   = zeros(Nx,Ny,Nz,3,3);

for g1=1:Nx
    for g2=1:Ny
        for g3=1:Nz
            PHM(g1,g2,g3,:,:)  = expm(squeeze(PHMa(g1,g2,g3,:,:)));
        end
    end
end

toc
PHM = ifftshift(PHM,1);
PHM = ifftshift(PHM,2);

%%
Uinx = zeros(Nx,Ny,Nview);
Uiny = zeros(Nx,Ny,Nview);
Uinz = zeros(Nx,Ny,Nview);

Ux   = zeros(Nx,Ny,Nview);
Uy   = zeros(Nx,Ny,Nview);
Uz   = zeros(Nx,Ny,Nview);
UIx  = zeros(Nx,Ny,Nview);
UIy  = zeros(Nx,Ny,Nview);
UIz  = zeros(Nx,Ny,Nview);

for Ang = 1:Nview
    Ang
    is_gauss        = 0;
    planar          = exp(1i*(kx_in(Ang)*X(:,:,1)*dx + ky_in(Ang)*Y(:,:,1)*dx + kz_in(Ang)*(Z(1,1,1)-1)*dz));  % Planar wave
    incident        = planar ./ exp(1i*k*(Z(1,1,1)-1)*dz);                              % Incident light

    if is_gauss == 1
        beam_fwhm_	= 50;                                                                     	% Full width half maximum [um]
        beam_scale_	= beam_fwhm_ / (2*sqrt(log(2)));                                            % Scale
        gaussian 	= exp(-(X(:,:,1)*dx/beam_scale_).^2 - (Y(:,:,1)*dx/beam_scale_).^2);        % Gaussian envelop
        incident   	= gaussian .* incident;    
    end

    incident        = ifftshift(incident);                             % Generating a Gaussian or planar beam, scalar
    incident        = ifft2(fft2(incident).*DFR_Ori);
    incident1       = 1/sqrt(2)*incident;                            % x
    incident2       = -1/sqrt(2)*incident;                           % y
    incident3       = zeros(Nx,Ny);                                    % z
    
    % Rotational Matrix Multiplication for input 4F system
    u_in1           = (Rin11(Ang)).*incident1+(Rin12(Ang)).*incident2+(Rin13(Ang)).*incident3;
    u_in2           = (Rin21(Ang)).*incident1+(Rin22(Ang)).*incident2+(Rin23(Ang)).*incident3;
    u_in3           = (Rin31(Ang)).*incident1+(Rin32(Ang)).*incident2+(Rin33(Ang)).*incident3;
    
    u_in_wo1        = u_in1;
    u_in_wo2        = u_in2;
    u_in_wo3        = u_in3;
    UU(:,:,Ang) = u_in_wo1;
    for g5=1:Nz    
        % Diffraction Part
        u_inx2 = ifft2(fft2(u_in1).*DFR.*(Numerical_Aperture));    
        u_iny2 = ifft2(fft2(u_in2).*DFR.*(Numerical_Aperture));
        u_inz2 = ifft2(fft2(u_in3).*DFR.*(Numerical_Aperture));
        % Phase Modulation Part
        u_inx3 = u_inx2.*PHM(:,:,g5,1,1)+u_iny2.*PHM(:,:,g5,1,2)+u_inz2.*PHM(:,:,g5,1,3);
        u_iny3 = u_inx2.*PHM(:,:,g5,2,1)+u_iny2.*PHM(:,:,g5,2,2)+u_inz2.*PHM(:,:,g5,2,3); 
        u_inz3 = u_inx2.*PHM(:,:,g5,3,1)+u_iny2.*PHM(:,:,g5,3,2)+u_inz2.*PHM(:,:,g5,3,3); 
        % Feedback Part
        u_in1  = u_inx3;
        u_in2  = u_iny3;
        u_in3  = u_inz3;
        
        % Without Sample
        u_in_wo1	= ifft2(fft2(u_in_wo1).*DFR);
        u_in_wo2	= ifft2(fft2(u_in_wo2).*DFR);
        u_in_wo3	= ifft2(fft2(u_in_wo3).*DFR);    
    end

    u_in1           = ifft2(fft2(u_in1).*DFR_Ori);
    u_in2           = ifft2(fft2(u_in2).*DFR_Ori);
    u_in3           = ifft2(fft2(u_in3).*DFR_Ori);
    u_in_wo1        = ifft2(fft2(u_in_wo1).*DFR_Ori);
    u_in_wo2        = ifft2(fft2(u_in_wo2).*DFR_Ori);
    u_in_wo3        = ifft2(fft2(u_in_wo3).*DFR_Ori);
    
    u_in1_F     = fft2(u_in1);
    u_in2_F     = fft2(u_in2);
    u_in3_F     = fft2(u_in3);
    u_in_wo1_F  = fft2(u_in_wo1);
    u_in_wo2_F  = fft2(u_in_wo2);
    u_in_wo3_F  = fft2(u_in_wo3);

    u1_F    = u_in1_F;
    u2_F    = u_in2_F;
    u3_F    = u_in3_F;
    ui1_F    = u_in_wo1_F;
    ui2_F    = u_in_wo2_F;
    ui3_F    = u_in_wo3_F;
    u1    = ifftshift(ifft2(u1_F));
    u2    = ifftshift(ifft2(u2_F));
    u3    = ifftshift(ifft2(u3_F));
    ui1   = ifftshift(ifft2(ui1_F));
    ui2   = ifftshift(ifft2(ui2_F));
    ui3   = ifftshift(ifft2(ui3_F));
    
    Ux(:,:,Ang) = u1;
    Uy(:,:,Ang) = u2;
    Uz(:,:,Ang) = u3;
    UIx(:,:,Ang)= ui1;
    UIy(:,:,Ang)= ui2;
    UIz(:,:,Ang)= ui3;

end

%%
Uxb = Ux;
Uyb = Uy;
Uzb = Uz;
UIxb = UIx;
UIyb = UIy;
UIzb = UIz;

%%
filename    = ['FW_BPM_',Sample_name, num2str(Nview), '_Ang_0_dn_', num2str(1000*wl), '_wl_', num2str(dx/wl), '_dx_Exp_Full3'];
save([Data_path, '\', filename, '.mat'], 'Uxa','Uya','Uza','Uxb','Uyb','Uzb','UIxa', 'UIya', 'UIza', 'UIxb', 'UIyb', 'UIzb', 'kx_in', 'ky_in', 'kz_in', 'dx', 'wl', 'n0', 'Nview', 'Nx', 'Ny', 'Nz', 'effective_NA', '-v7.3');
disp(filename);
