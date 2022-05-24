%% 05/25/12
%% Seyedalireza Abootorabi
%% This code uses the Reynolds stress to compute the velocity profile unlike the other version which uses the eddy-viscosity
%%=========================================================================
%% model-based control of turbulent flow
%%
%% channel flow 
%%
%% transverse wall oscillations on both walls (in phase)
%% W(y = \pm 1,t) = 2*epsilon*sin(wt t)
%%
%% NS equations linearized around the turbulent mean velocity that is 
%% obtained under the assumption that the turbulent viscosity if the flow
%% with no control captures the influence of background turbulence on the
%% mean velocity
%%
%% determine correction to drag reduction (saved power) directly from
%% the second-order correction to "uv".
%%=========================================================================

%%=========================================================================
%% SUMMARY:
%%=========================================================================
%% SECTION 1: 
%% load the turbulent mean velocity and the Reynolds-stress 
%% obtained from DNS of del Alamo and Jimenez JFM 2003 and 2004.
%%=========================================================================
%% SECTION 2: 
%% set input parameters
%%=========================================================================
%% SECTION 3: 
%% initialize the code
%%=========================================================================
%% SECTION 4: 
%% Cess turbulent eddy viscosity
%% Reynolds-Tiederman turbulent mean velocity in streamwise direction
%%=========================================================================
%% SECTION 5: 
%% extrapolate DNS data
%%=========================================================================
%% SECTION 6: 
%% compute the turbulent energy spectrum, \bar{E} (kappa), obtained from DNS
%%=========================================================================
%% SECTION 7:
%% spanwise turbulent mean velocity, W0 (y,t), obtained from nuT of the 
%% flow with no control
%%=========================================================================
%% SECTION 8:
%% required power, Preq0, for maintaining the transverse oscillations under
%% the assumption that nuT of the flow with no control determines the
%% effect of background turbulence
%%=========================================================================
%% SECTION 9:
%% linearized equations around base flow: [U0(y), 0, W0(y,t)]
%%=========================================================================
%% SECTION 10:
%% obtain the second-order correction to uv using Lyapunov equations
%%=========================================================================
%% SECTION 11:
%% obtain the second-order corrections to turbulent mean velocity and the
%% saved power
%%=========================================================================
%% SECTION 12:
%% integrate over horizontal wavenumbers to obtain total
%% effect of control on each quantity 
%% NOTE: in the initial submission of our JFM 2011 paper, we have reported 
%% the integrated results
%%=========================================================================
%% SECTION 13:
%% plot results
%%=========================================================================

clear
format compact
format short e
clc

%%=========================================================================
%% SECTION 1: 
%% load the turbulent mean velocity and the Reynolds-stress 
%% obtained from DNS of del Alamo and Jimenez JFM 2003 and 2004.
%%=========================================================================

% load turbulent mean profiles (averaged over x and z)
load DelAlamo_Jimenez_JFM2004_mean_profiles_R186
% load('/home/zarex004/Documents/MATLAB/Turbulent_flow/DelAlamo_Jimenez_JFM2004_Data/DelAlamo_Jimenez_JFM2004_mean_profiles_R186')
% load DelAlamo_Jimenez_JFM2004_mean_profiles_R186

y_DNS  = y;
U0_DNS = Umean;
uu_DNS = u_rms.^2;
vv_DNS = v_rms.^2;
ww_DNS = w_rms.^2;
uv_DNS = uv;
clear Umean omegax_rms omegay_rms omegaz_rms u_rms v_rms w_rms uv uw vw y yp

% load turbulent spectrum at different horizontal wavenumbers and a limited
% number of points in wall-normal direction
load DelAlamo_Jimenez_JFM2004_spectra_uvw_R186
% load('/home/zarex004/Documents/MATLAB/Turbulent_flow/DelAlamo_Jimenez_JFM2004_Data/DelAlamo_Jimenez_JFM2004_spectra_uvw_R186')
% load DelAlamo_Jimenez_JFM2004_spectra_uvw_R186

% set mean mode (kx = 0, kz = 0) to zero
uu_spectra_DNS(:,1,1) = 0;
vv_spectra_DNS(:,1,1) = 0;
ww_spectra_DNS(:,1,1) = 0;
uv_spectra_DNS(:,1,1) = 0;

%%=========================================================================
%% SECTION 2: 
%% set input parameters
%%=========================================================================

% friction Reynolds number

R = 186;

% number of collocation points inwall-normal direction

N = 31;

% the horizontal wavenumbers bkappa = (kx,kz)
% these values have to be within the range of values that del Alamo and Jimenez have selected
% for the interpolation to give correct results, set the lower and upper
% bounds to be, respectively, just above or below the range of values used
% by del Alamo and Jimenez. Otherwise, interpolation gives NaN at the
% boundary of the interval.

% wavenumber in the streamwise direction
kxval = logspace(log10(kxval_DNS(2)),log10(kxval_DNS(end)-0.01),50)';
d1 = log(kxval(2)/kxval(1));
kxval1 = kxval(1)/exp(d1)^25;
kxvalnew = logspace(log10(kxval1),log10(kxval_DNS(end)-0.01),75)';

% kxval = kxvalnew(1:25);
kxval = kxvalnew;
kxgrd  = length(kxval);

% wavenumber in the spanwise driection
kzval = logspace(log10(kzval_DNS(2)),log10(kzval_DNS(end)-0.01),51)';
kzgrd = length(kzval);

% the period of oscillations in viscous units, T^+

%Tplusval = linspace(30,500,5);
Tplusval = [30:60:300];
% Tplusval = 102.5;

%%=========================================================================
%% SECTION 3: 
%% initialize the code
%%=========================================================================

% the frequency of oscillations in outer units

wtval = 2*pi*R./Tplusval;
wtgrd = length(wtval);

% the unit imaginary number

ii = sqrt(-1);

% differentiation matrices

[yT,DM] = chebdif(N+2,3);

D1T = DM(:,:,1);
D2T = DM(:,:,2);
D3T = DM(:,:,3);

% implement homogeneous boundary conditions

D1 = DM(2:N+1,2:N+1,1);
D2 = DM(2:N+1,2:N+1,2);
D3 = DM(2:N+1,2:N+1,3);

% fourth derivative with clamped conditions

[y,D4]=cheb4c(N+2);

% integration weights

IWT = INTweights(N+2,2);
IW = IWT(2:end-1,2:end-1);
invIW = inv(IW);

% identity and zero matrices

I = eye(N);
Z = zeros(N,N);
Ibig = eye(2*N);

% initialize variables by zeros

uv_2 = zeros(N,kxgrd,kzgrd,wtgrd);

uvvec_2 = zeros(N,wtgrd);

%%=========================================================================
%% SECTION 4: 
%% Cess turbulent eddy viscosity
%% Reynolds-Tiederman turbulent mean velocity in streamwise direction
%%=========================================================================

% set the parameters in the Cess turbulent eddy viscosity

if R == 186
    kappa = 0.61;  
    Amp   = 46.2;  
elseif R == 547
    kappa = 0.455; 
    Amp   = 29.4;  
elseif R == 934
    kappa = 0.43;  
    Amp   = 27;
else
    display('Using nuT with parameters fitted for R = 2000!')
    kappa = 0.426; 
    Amp   = 25.4;  
end

% Cess turbulent viscosity profile nuT

nuT = 0.5*( (1 + ( (1/3)*kappa*R*(1 - y.^2).*(1 + 2*y.^2).*(1 - exp(-(1 - abs(y))*R/Amp)) ).^2 ).^(1/2) - 1  );
nuTT = 0.5*( (1 + ( (1/3)*kappa*R*(1 - yT.^2).*(1 + 2*yT.^2).*(1 - exp(-(1 - abs(yT))*R/Amp)) ).^2 ).^(1/2) - 1  );

% derivatives of nuT

nuTyT = D1T*nuTT;
nuTyyT = D2T*nuTT;
nuTyyyT = D3T*nuTT;

nuTy = nuTyT(2:end-1);
nuTyy = nuTyyT(2:end-1);
nuTyyy = nuTyyyT(2:end-1);

% base flow [U(y), 0, W(y,t)]
% streamwise turbulent mean velocity in the flow with no control U0 (y)

LHS = diag(1+nuT)*D2 + diag(nuTy)*D1;
invLHS = inv(LHS);

U0 = LHS\(-R*ones(N,1));
U0y = D1*U0;
U0yy = D2*U0;

% the bulk flux in the flow with no control

UB = sum(IW'*U0)/2;

%%=========================================================================
%% SECTION 5: 
%% extrapolate DNS data 
%%=========================================================================

% extrapolate mean profiles on the wall-normal collocation points

uu_DNS = interp1(y_DNS,uu_DNS,y);
vv_DNS = interp1(y_DNS,vv_DNS,y);
ww_DNS = interp1(y_DNS,ww_DNS,y);
uv_DNS = interp1(y_DNS,uv_DNS,y);


% extrapolate spectra on the wall-normal collocation points and on the horizontal wavenumbers

% find the closest working y to the y_spectra used in DNS

for indy = 1:Ny_spectra
    [m indy_spectra(indy,1)] = min(abs(y-y_spectra_DNS(indy)));
end

% add the wall (y == 0) to the DNS data

uu_spectra_DNS_temp = zeros(Ny_spectra+1,length(kxval_DNS),length(kzval_DNS));
vv_spectra_DNS_temp = zeros(Ny_spectra+1,length(kxval_DNS),length(kzval_DNS));
ww_spectra_DNS_temp = zeros(Ny_spectra+1,length(kxval_DNS),length(kzval_DNS));
uv_spectra_DNS_temp = zeros(Ny_spectra+1,length(kxval_DNS),length(kzval_DNS));
uu_spectra_DNS_temp(2:end,:,:) = uu_spectra_DNS;
vv_spectra_DNS_temp(2:end,:,:) = vv_spectra_DNS;
ww_spectra_DNS_temp(2:end,:,:) = ww_spectra_DNS;
uv_spectra_DNS_temp(2:end,:,:) = uv_spectra_DNS;

% extrapolate

[yi_spectra_DNS,kxvali_DNS,kzvali_DNS] = ndgrid([-1; y_spectra_DNS],kxval_DNS,kzval_DNS);
[yi,kxvali,kzvali] = ndgrid(y(indy_spectra),kxval,kzval);
uu_spectra_DNS = interpn(yi_spectra_DNS,kxvali_DNS,kzvali_DNS,uu_spectra_DNS_temp,yi,kxvali,kzvali,'linear',0);
vv_spectra_DNS = interpn(yi_spectra_DNS,kxvali_DNS,kzvali_DNS,vv_spectra_DNS_temp,yi,kxvali,kzvali,'linear',0);
ww_spectra_DNS = interpn(yi_spectra_DNS,kxvali_DNS,kzvali_DNS,ww_spectra_DNS_temp,yi,kxvali,kzvali,'linear',0);
uv_spectra_DNS = interpn(yi_spectra_DNS,kxvali_DNS,kzvali_DNS,uv_spectra_DNS_temp,yi,kxvali,kzvali,'linear',0);
clear uu_spectra_DNS_temp vv_spectra_DNS_temp ww_spectra_DNS_temp uv_spectra_DNS_temp

%%=========================================================================
%% SECTION 6: 
%% compute the turbulent energy spectrum, \bar{E} (bkappa), obtained from DNS
%%=========================================================================

for indkx = 1:kxgrd
    for indkz = 1:kzgrd 
        
        % interpolate DNS-based spectrum on y
        
        uu_temp = uu_spectra_DNS(:,indkx,indkz);
        uu_temp = interp1([-1; y_spectra_DNS],[0; uu_temp],y(floor(N/2)+1:end),'linear',0);
        vv_temp = vv_spectra_DNS(:,indkx,indkz);
        vv_temp = interp1([-1; y_spectra_DNS],[0; vv_temp],y(floor(N/2)+1:end),'linear',0);
        ww_temp = ww_spectra_DNS(:,indkx,indkz);
        ww_temp = interp1([-1; y_spectra_DNS],[0; ww_temp],y(floor(N/2)+1:end),'linear',0);
        
        % E (y,bkappa), on the lower half of the channel
        
        E_turb_profile(:,indkx,indkz) = uu_temp + vv_temp + ww_temp;
        
        % integrate in y to get energy spectrum \bar{E} (bkappa)
        % integrate over the lower half and then multiply by 2
        
        E_turb(indkx,indkz) = 2*diag(IW(floor(N/2)+1:end,floor(N/2)+1:end))'*(E_turb_profile(:,indkx,indkz));
    end
end
clear uu_temp vv_temp ww_temp

% E (y,bkappa) on both lower and upper halves

E_turb_profile = [E_turb_profile(end:-1:2,:,:); E_turb_profile(:,:,:)];


% %% load noise-modelling covariances
% load('/Users/arminzare/Documents/MATLAB/LNSE/NoiseModeling/LowRank_Fu/Results/NoiseModelling/Qmat_data_N51_75*51.mat')
% % load('/home/zarex004/Documents/MATLAB/NoiseModeling/LowRank_Fu/Results/NoiseModelling/Qmat_data_N51_75*51.mat')

% loop over different values for frequency of oscillations

for indwt = 1:wtgrd
    
    wt = wtval(indwt);

    %%=====================================================================
    %% SECTION 7:
    %% spanwise turbulent mean velocity, W0 (y,t), obtained from nuT of the 
    %% flow with no control
    %%=====================================================================
    
    % compute W0 (y,t)
    % W0 (y,t) = alpha (Wp (y) e^{ii*w*t} + Wm (y) e^{-ii*w*t})
    
    % Wp
    
    Wp1bar = (LHS/R-ii*wt*I)\(wt*ones(N,1));
    Wp1ybar = D1*Wp1bar;
    Wp1yybar = D2*Wp1bar;
    Wp1(:,indwt) = Wp1bar + (-ii)*ones(N,1);
    Wp1y(:,indwt) = Wp1ybar;
    Wp1yy(:,indwt) = Wp1yybar;

    % Wm
    
    Wm1bar = (LHS/R+ii*wt*I)\(wt*ones(N,1));
    Wm1ybar = D1*Wm1bar;
    Wm1yybar = D2*Wm1bar;
    Wm1(:,indwt) = Wm1bar + (ii)*ones(N,1);
    Wm1y(:,indwt) = Wm1ybar;
    Wm1yy(:,indwt) = Wm1yybar;

    %%=====================================================================
    %% SECTION 8:
    %% required power, Preq0, for maintaining the transverse oscillations 
    %% under the assumption that nuT of the flow with no control determines 
    %% the effect of background turbulence
    %%=====================================================================

    Wp1Tybar = D1T*[0; Wp1bar; 0];
    Wp1Ty = Wp1Tybar;
    Wm1Tybar = D1T*[0; Wm1bar; 0];
    Wm1Ty = Wm1Tybar;
    Preq_0(indwt) = 100/(2*UB)*(imag(Wp1Ty(end) - Wm1Ty(end))/R - ...
                                imag(Wp1Ty(1) - Wm1Ty(1))/R);
                
    % loop over streamwise wavenumber
    
    for indkx = 1:kxgrd
        
        kx = kxval(indkx);
        
        % loop over spanwise wavenumber
    
        for indkz = 1:kzgrd

            kz = kzval(indkz);

            k2 = kx*kx + kz*kz;
            k4 = k2*k2;

            %%=============================================================
            %% SECTION 9:
            %% linearized equations around base flow: [U0(y), 0, W0(y,t)]
            %% [v eta] state variables
            %%=============================================================

            % Laplacian
            
            Delta = (D2 - k2*I);
            Delta2 = (D4 - 2*k2*D2 + k4*I);

            % integration weights that relate the iner-product in Hos+L2 to L2
            QvT = (1/k2)*D1T'*IWT*D1T + IWT;
            QetaT = (1/k2)*IWT;
            
            Qv = QvT(2:N+1,2:N+1);
            Qeta = QetaT(2:N+1,2:N+1);
            
            % square roots of Qv and Qeta
            Qv_half = sqrtm(Qv);
            Qeta_half = sqrtm(Qeta);
            
            Qhalf = [Qv_half, Z; Z, Qeta_half];
            invQhalf = [Qv_half\I, Z; Z, Qeta_half\I];
            
            % inverse of Q
            invQv = Qv\I;
            invQeta = Qeta\I;
            Q = [Qv, Z; Z, Qeta];
            invQ = [invQv, Z; Z, invQeta];
            
            
            % A0 operator (order of alpha^0)
            A11 = Delta\(1/R*(diag(1+nuT))*Delta2 + 1/R*diag(nuTyy)*(D2+k2*I) + 2/R*diag(nuTy)*(D3 - k2*D1) + ...
                  ii*kx*(diag(U0yy) - diag(U0)*Delta));
            A12 = Z;
            A21 = -ii*kz*diag(U0y);
            A22 = 1/R*diag(1+nuT)*Delta + 1/R*diag(nuTy)*D1 - ii*kx*diag(U0);

            A0 = [A11, A12; A21, A22];
            
            % A0 in new coordiantes
            A0 = Qhalf*A0*invQhalf;
            
            % Am1 operator (order of alpha^1)
            A11 = Delta\(ii*kz*(diag(Wm1yy(:,indwt)) - diag(Wm1(:,indwt))*Delta));
            A12 = Z;
            A21 = ii*kx*diag(Wm1y(:,indwt));
            A22 = -ii*kz*diag(Wm1(:,indwt));

            Am1 = [A11, A12; A21, A22];
            
            % Am1 in new coordiantes
            Am1 = Qhalf*Am1*invQhalf;
            
            % Ap1 operator (order of alpha^1)
            A11 = Delta\(ii*kz*(diag(Wp1yy(:,indwt)) - diag(Wp1(:,indwt))*Delta));
            A12 = Z;
            A21 = ii*kx*diag(Wp1y(:,indwt));
            A22 = -ii*kz*diag(Wp1(:,indwt));

            Ap1 = [A11, A12; A21, A22];
            
            % Ap1 in new coordiantes
            Ap1 = Qhalf*Ap1*invQhalf;
            
            %%=============================================================
            %% SECTION 10:
            %% obtain the second-order statistics of velocity using 
            %% Lyapunov equations
            %%=============================================================

            % P0 (order of alpha^0)
            
            % force with spectrum of the turbulent flow
            % based on intuition frm Homogeneous Isotropic Turbulence
            % forcing: Eturb^(1/2) (Eturb^{ad})^(1/2) = Eturb^(1/2) * Q^{-1} * Eturb^(1/2) * Q
            % its kernel representation is Eturb^(1/2) * Q^{-1} * Eturb^(1/2)
            
            F0 = A0;
            
            P0 = lyap(F0,[diag(E_turb_profile(:,indkx,indkz)), Z; ...
                             Z, diag(E_turb_profile(:,indkx,indkz))]);
            
            % energy in the flow with no control             
            
            E_0(indkx,indkz,indwt) = trace(P0);
            
            % multiply Pbar0 by a scalar such that the energy spectrum
            % obtained from the model is equal to the energy obtained from
            % DNS. i.e. myscalar = E_turb/E_0
            
            myscalar = E_turb(indkx,indkz)/E_0(indkx,indkz,indwt);
            sigma_shape(indkx,indkz) = myscalar;
            
            P0 = P0*myscalar;
            E_0(indkx,indkz,indwt) = E_0(indkx,indkz,indwt)*myscalar;
            
            % P1 (order of alpha^1)

            Fp1 = A0 + ii*wt*Ibig;
            
            RHS = Am1*P0 + ...
                  P0*Ap1';
            
            P1 = lyap(Fp1,F0',RHS);

            % P2 (order of alpha^2)
            
            RHS = Ap1*P1 + Am1*P1' + ...
                  P1'*Ap1' + P1*Am1';

            P20 = lyap(F0,RHS);

            % second-oredr correction in alpha to the energy in the flow with control
            
            E_2(indkx,indkz,indwt) = trace(P20);

            %%=============================================================
            %% SECTION 10:
            %% obtain the second-order correction to uv using Lyapunov equations
            %%=============================================================

            % the C operator mapping state to velocity field
            
            Cu = [ii*kx*D1, -ii*kz*I]/k2;
            Cv = [I, Z];
            Cw = [ii*kz*D1, ii*kx*I]/k2;
            
            % 2nd-order correction to uv in the flow with control
            
            uv_2(:,indkx,indkz,indwt) = real(diag(Cu*invQhalf*P20*invQhalf*Cv'));
            
            %%=============================================================
            %% SECTION 11:
            %% obtain the second-order corrections to turbulent mean 
            %% velocity and the saved power
            %%=============================================================
            
            % the 2nd-order correction to U if the pressure gradient is
            % kept constant and the bulk flux is allowed to change
            
            U_2(:,indkx,indkz,indwt) = D2\(R*(D1*uv_2(:,indkx,indkz,indwt)));

            % the 2nd-order correction to the bulk flux if the pressure 
            % gradient is kept constant and the bulk flux is allowed to change
            
            UB_2(indkx,indkz,indwt) = diag(IW)'*U_2(:,indkx,indkz,indwt)/2;
            
            % the 2nd-order correction to pressure gradient for constant bulk flux
            % in our JFM 2011 paper, this is denoted by -\tau_{w2}
            % this is equal to the saved power
            % Psave = alpha^2 Px_nuT2
            
            Px_2(indkx,indkz,indwt) = UB_2(indkx,indkz,indwt)/(UB);
            
            % the 2nd-order correction to U if the pressure gradient is 
            % adjusted to obtain constant bulk flux
            
            U_2_UBconst(:,indkx,indkz,indwt) = U_2(:,indkx,indkz,indwt) - U0*Px_2(indkx,indkz,indwt);
            
            [indwt indkx indkz]

        end
    end
end

rrr

%%=========================================================================
%% SECTION 12:
%% integrate over horizontal wavenumbers to obtain total
%% effect of control on each quantity 
%% NOTE: in the initial submission of our JFM 2011 paper, we have reported 
%% the integrated results
%%=========================================================================

% integration scale for integrating the functions in log scale in kx and kz

logkxkz_scale = log(kxval(2)/kxval(1))*log(kzval(2)/kzval(1));

% premultiply by kx*kz for integration in log scale

for indkx = 1:kxgrd
    
    kx = kxval(indkx);
    
    for indkz = 1:kzgrd
        
        kz = kzval(indkz);
        
        uv_2_pre(:,indkx,indkz,:)        = kx*kz*uv_2(:,indkx,indkz,:);
        U_2_pre(:,indkx,indkz,:)         = kx*kz*U_2(:,indkx,indkz,:);
        U_2_UBconst_pre(:,indkx,indkz,:) = kx*kz*U_2_UBconst(:,indkx,indkz,:);
        UB_2_pre(indkx,indkz,:) = kx*kz*UB_2(indkx,indkz,:);
        Px_2_pre(indkx,indkz,:) = kx*kz*Px_2(indkx,indkz,:);
        
        E_0_pre(indkx,indkz,:) = kx*kz*E_0(indkx,indkz,:);
        E_2_pre(indkx,indkz,:) = kx*kz*E_2(indkx,indkz,:);
    end
    
end

            
% integrate over horizontal wavenumbers

for indwt = 1:wtgrd
    
    uv_2_avg(:,indwt)        = sum(sum(uv_2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    U_2_avg(:,indwt)         = sum(sum(U_2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    U_2_UBconst_avg(:,indwt) = sum(sum(U_2_UBconst_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    UB_2_avg(indwt) = sum(sum(UB_2_pre(:,:,indwt),2),1)*logkxkz_scale;
    Px_2_avg(indwt) = sum(sum(Px_2_pre(:,:,indwt),2),1)*logkxkz_scale;
    
end

clear kxvali_DNS kzvali_DNS yi_spectra_DNS

cd Results
save Results_TWO_pert_uv_Rashad_newcoordinates_N51_Tfull

cd ..

error('This error is just to terminate the code. Please proceed to plotting the results! - Rashad')

%%=========================================================================
%% SECTION 13:
%% plot results
%%=========================================================================

% interpolate data for high-quality figures
% in our paper, I used enough points in T^+ to capture the trends, interpolation is used only to obtain a smoother curve

% use T^+ (the 0.01 is used to avoid NaN after interpolation)
Tplus_plot = (Tplusval(1)+0.01):0.5:(Tplusval(end)-0.01);


% uv
figure(11); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,uv_2_NM1_avg,yvec/R-1,'pchip'),'r-');

% Psave2
figure(1); hold on; plot(Tplus_plot,100*interp1(Tplusval,Px_2_avg,Tplus_plot,'spline'),'k-');
figure(2); hold on; f = interp1(Tplusval,Px_2_avg,Tplus_plot,'spline'); plot(Tplus_plot,f/max(f),'k-');

% compare with DNS data: Quadrio and Ricco, JFM 2004

load('/Users/rashad/Backups/Rashad Office Desktop/Rashad/Research data/Quadrio-Ricco-jfm2004-data/Quadrio-Ricco-jfm2004-data.mat')
% Wp_4p5, Wp_12, and Wp_18 are DNS data from Quadrio and Ricco, JFM 2004
% example: 
% Wp_4p5 is the data for alpha = 4.5 in their paper (alpha = 4.5/2 = 2.25 in our paper)
% Wp_4p5 has three columns:
% first column: the values of T^+
% second column: saved power
% third column: required power 
% note: they define a negative required power, we defined positive required power
% for them: Pnet = Psave + Preq
% for us:   Pnet = Psave - Preq

figure(2); hold on; f = Wp_4p5(:,2); plot(Wp_4p5(:,1),f/max(f),'ko');
figure(2); hold on; f = Wp_12(:,2); plot(Wp_12(:,1),f/max(f),'ks');
figure(2); hold on; f = Wp_18(:,2); plot(Wp_18(:,1),f/max(f),'kv');

% Preq0
figure(3); hold on; plot(Tplus_plot,interp1(Tplusval,Preq_0,Tplus_plot,'spline'),'k-');
% compare with DNS data: Quadrio and Ricco, JFM 2004
figure(4); hold on; f = interp1(Tplusval,Preq_0,Tplus_plot,'spline'); plot(Tplus_plot,f/max(f),'k-');
figure(4); hold on; f = Wp_4p5(:,3); plot(Wp_4p5(:,1),f/f(2),'ko');
figure(4); hold on; f = Wp_12(:,3); plot(Wp_12(:,1),f/f(1),'ks');
figure(4); hold on; f = Wp_18(:,3); plot(Wp_18(:,1),f/f(1),'kv');

% Pnet2
figure(7); hold on; plot(Tplus_plot,100*interp1(Tplusval,Px_2_avg,Tplus_plot,'spline')-interp1(Tplusval,Preq_0,Tplus_plot,'spline'),'k-');
load('/Users/rashad/Backups/Rashad Office Desktop/Rashad/Research data/Quadrio-Ricco-jfm2004-data/Quadrio-Ricco-jfm2004-data.mat')
figure(8); hold on; f = 100*interp1(Tplusval,Px_2_avg,Tplus_plot,'spline')-interp1(Tplusval,Preq_0,Tplus_plot,'spline'); plot(Tplus_plot,f/max(f),'k-');
figure(8); hold on; f = Wp_4p5(:,2)+Wp_4p5(:,3); plot(Wp_4p5(:,1),f/max(f),'ko');
figure(18); hold on; f = 100*interp1(Tplusval,Px_2_avg,Tplus_plot,'spline')-interp1(Tplusval,Preq_0,Tplus_plot,'spline'); plot(Tplus_plot,f,'k-');
figure(18); hold on; f = Wp_4p5(:,2)+Wp_4p5(:,3); plot(Wp_4p5(:,1),f/2.25^2,'ko');

% U0 and nuT0
figure(9); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,U0,yvec/R-1,'pchip'),'k-');
figure(10); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,nuT,yvec/R-1,'pchip'),'k-');

% U2
figure(20); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,U_2_UBconst_avg,yvec/R-1,'spline'),'k-');

% Wp0
figure(22); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,real(Wp1),yvec/R-1,'spline'),'k-');
figure(23); hold on; yvec = logspace(log10(0.1),log10(R),100); plot(yvec,interp1(y,imag(Wp1),yvec/R-1,'spline'),'k-');
