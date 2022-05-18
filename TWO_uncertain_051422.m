%% 05/13/2022
%%=========================================================================
% 
% Written by: Arya Abootorabi and Armin Zare
% TWO code from moajov'12 in the presence of parametric unertainty
%
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
%% U2date turbulent viscosity using nu_T = c k^2/epsilon
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
%% obtain the second-order statistics of velocity using Lyapunov equations
%%=========================================================================
%% SECTION 11:
%% obtain the second-order corrections to k, epsilon, and \nuT
%%=========================================================================
%% SECTION 12:
%% obtain the second-order corrections to turbulent mean velocity and the
%% required power and the saved power
%%=========================================================================
%% SECTION 13:
%% integrate k, epsilon, etc over horizontal wavenumbers to obtain total
%% effect of control on each quantity 
%% NOTE: in the initial submission of our JFM 2011 paper, we have reported 
%% the integrated results
%%=========================================================================
%% SECTION 14:
%% plot results
%%=========================================================================

clear
clc
format compact

%%=========================================================================
%% SECTION 1: 
%% load the turbulent mean velocity and the Reynolds-stress 
%% obtained from DNS of del Alamo and Jimenez JFM 2003 and 2004.
%%=========================================================================

% load turbulent mean profiles (averaged over x and z)

load DelAlamo_Jimenez_JFM2004_mean_profiles_R186.mat
y_DNS  = y;
U0_DNS = Umean;
uu_DNS = u_rms.^2;
vv_DNS = v_rms.^2;
ww_DNS = w_rms.^2;
uv_DNS = uv;
clear Umean omegax_rms omegay_rms omegaz_rms u_rms v_rms w_rms uv uw vw y yp

% load turbulent spectrum at different horizontal wavenumbers and a limited
% number of points in wall-normal direction

load DelAlamo_Jimenez_JFM2004_spectra_uvw_R186.mat

%%=========================================================================
%% SECTION 2: 
%% set input parameters
%%=========================================================================

% friction Reynolds number

R = 186;

% number of collocation points inwall-normal direction

N = 51;

% wave amplitude
alpha = 2.25;

% imaginary unit
ii = sqrt(-1);

% identity and zero matrices
I = eye(N);
Z = zeros(N,N);
I2 = eye(2*N); 
ZZ = zeros(2*N, 2*N);

% parameters for uncertainty gamma_alpha, gamma_theta
% with uncertainty
sigma2_alpha = 30;
%sigma2_theta = 0.69;
sigma2_theta = 0.005;

mu_alpha = 0;
mu_theta = 0;

mup1 = (1 + mu_alpha)*exp(ii*mu_theta + sigma2_theta/2) - 1;
mum1 = (1 + mu_alpha)*exp(-ii*mu_theta + sigma2_theta/2) - 1;

sigma2p2 = (exp(2*ii*mu_theta + sigma2_theta))*(sigma2_alpha*exp(sigma2_theta) + ((1 + mu_alpha)^2)*(exp(sigma2_theta) - 1));
sigma2m2 = (exp(-2*ii*mu_theta + sigma2_theta))*(sigma2_alpha*exp(sigma2_theta) + ((1 + mu_alpha)^2)*(exp(sigma2_theta) - 1));

            
            
% the horizontal wavenumbers bkappa = (kx,kz)
% these values have to be within the range of values that del Alamo and Jimenez have selected
% for the interpolation to give correct results, set the lower and U2per
% bounds to be, respectively, just above or below the range of values used
% by del Alamo and Jimenez. Otherwise, interpolation gives NaN at the
% boundary of the interval.

% wavenumber in the streamwise direction

kxval = logspace(log10(kxval_DNS(1)+0.01),log10(kxval_DNS(end)-0.01),50)';
kxgrd  = length(kxval);

% wavenumber in the spanwise driection

kzval = logspace(log10(kzval_DNS(1)+0.01),log10(kzval_DNS(end)-0.01),51)';
kzgrd  = length(kzval);

% the period of oscillations in viscous units, T^+

% Tplusval = linspace(10,300,10);
Tplusval = 102.5;

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

[yvecT,DM] = chebdif(N+2,3);
yvec=yvecT(2:end-1);

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
IW2 = kron(eye(2),IW);
IW3 = kron(eye(3),IW);
IWhalf = sqrtm(IW);
invIWhalf = IWhalf\I;
IW3half = sqrtm(IW3);
invIW3half = IW3half\IW3;

% streamwise mean velocity
U = diag(1-yvec.^2);
Uy = diag(-2*yvec);
Uyy = diag(-2*ones(size(yvec)));
% identity and zero matrices

I = eye(N);
Z = zeros(N,N);
Ibig = eye(2*N);

% initialize variables by zeros

uu_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uv_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uw_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vv_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vw_0 = zeros(N,kxgrd,kzgrd,wtgrd);
ww_0 = zeros(N,kxgrd,kzgrd,wtgrd);

uu_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uv_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uw_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vv_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vw_2 = zeros(N,kxgrd,kzgrd,wtgrd);
ww_2 = zeros(N,kxgrd,kzgrd,wtgrd);

uuvec_0 = zeros(N,wtgrd);
uvvec_0 = zeros(N,wtgrd);
uwvec_0 = zeros(N,wtgrd);
vvvec_0 = zeros(N,wtgrd);
vwvec_0 = zeros(N,wtgrd);
wwvec_0 = zeros(N,wtgrd);

uuvec_2 = zeros(N,wtgrd);
uvvec_2 = zeros(N,wtgrd);
uwvec_2 = zeros(N,wtgrd);
vvvec_2 = zeros(N,wtgrd);
vwvec_2 = zeros(N,wtgrd);
wwvec_2 = zeros(N,wtgrd);

uxux_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vyvy_0 = zeros(N,kxgrd,kzgrd,wtgrd);
wzwz_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uyvx_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uzwx_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vzwy_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uyuy_0 = zeros(N,kxgrd,kzgrd,wtgrd);
wywy_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vxvx_0 = zeros(N,kxgrd,kzgrd,wtgrd);
wxwx_0 = zeros(N,kxgrd,kzgrd,wtgrd);
uzuz_0 = zeros(N,kxgrd,kzgrd,wtgrd);
vzvz_0 = zeros(N,kxgrd,kzgrd,wtgrd);

uxux_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vyvy_2 = zeros(N,kxgrd,kzgrd,wtgrd);
wzwz_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uyvx_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uzwx_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vzwy_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uyuy_2 = zeros(N,kxgrd,kzgrd,wtgrd);
wywy_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vxvx_2 = zeros(N,kxgrd,kzgrd,wtgrd);
wxwx_2 = zeros(N,kxgrd,kzgrd,wtgrd);
uzuz_2 = zeros(N,kxgrd,kzgrd,wtgrd);
vzvz_2 = zeros(N,kxgrd,kzgrd,wtgrd);


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

% E (y,bkappa) on both lower and U2per halves

E_turb_profile = [E_turb_profile(end:-1:2,:,:); E_turb_profile(:,:,:)];

% nuT, k, and epsilon in the flow with no control

% nuT is the Cess profile

nuT_DNS = nuT;

% k is obtained from DNS

k_DNS = (uu_DNS + vv_DNS + ww_DNS)/2;
k_0RT = k_DNS;

% epsilon is obtained from the k-epsilon model

epsilon_DNS = 0.09*R^2*k_DNS.^2./nuT_DNS;
epsilon_0RT = epsilon_DNS;


% pre-allocations
E_0 = zeros(kxgrd,kzgrd,wtgrd);
E_2 = zeros(kxgrd,kzgrd,wtgrd);
pre_E2 = zeros(kxgrd,kzgrd,wtgrd);
E = zeros(kxgrd,kzgrd,wtgrd);

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
            
            QvT = (D1T'*IWT*D1T/k2 + IWT);
            Qv = QvT(2:N+1,2:N+1);
            invQv = inv(Qv);
            Qeta = IW/k2;
            invQeta = inv(Qeta);

            Q = [Qv, Z; Z, Qeta];
            invQ = [invQv, Z; Z, invQeta];
            
            % A0 operator (order of alpha^0)
            
            A11 = Delta\(1/R*(diag(1+nuT))*Delta2 + 1/R*diag(nuTyy)*(D2+k2*I) + 2/R*diag(nuTy)*(D3 - k2*D1) + ...
                  ii*kx*(diag(U0yy) - diag(U0)*Delta));
            A12 = Z;
            A21 = -ii*kz*diag(U0y);
            A22 = 1/R*diag(1+nuT)*Delta + 1/R*diag(nuTy)*D1 - ii*kx*diag(U0);

            A0 = [A11, A12; A21, A22];

            % Am1 operator (order of alpha^1)
            
            A11 = Delta\(ii*kz*(diag(Wm1yy(:,indwt)) - diag(Wm1(:,indwt))*Delta));
            A12 = Z;
            A21 = ii*kx*diag(Wm1y(:,indwt));
            A22 = -ii*kz*diag(Wm1(:,indwt));

            Am1 = [A11, A12; A21, A22];

            % Ap1 operator (order of alpha^1)
            
            A11 = Delta\(ii*kz*(diag(Wp1yy(:,indwt)) - diag(Wp1(:,indwt))*Delta));
            A12 = Z;
            A21 = ii*kx*diag(Wp1y(:,indwt));
            A22 = -ii*kz*diag(Wp1(:,indwt));

            Ap1 = [A11, A12; A21, A22];

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
            
            Pbar0 = lyap(F0,sqrt([diag(E_turb_profile(:,indkx,indkz)), Z; ...
                             Z, diag(E_turb_profile(:,indkx,indkz))])* ...
                             invQ* ...
                             sqrt([diag(E_turb_profile(:,indkx,indkz)), Z; ...
                             Z, diag(E_turb_profile(:,indkx,indkz))]));
                        
            
            % energy in the flow with no control             
            
            E_0(indkx,indkz,indwt) = trace(Pbar0*Q);
            
            % multiply Pbar0 by a scalar such that the energy spectrum
            % obtained from the model is equal to the energy obtained from
            % DNS. i.e. myscalar = E_turb/E_0
            
            myscalar = E_turb(indkx,indkz)/E_0(indkx,indkz,indwt);
            sigma_shape(indkx,indkz) = myscalar;
            
            Pbar0 = Pbar0*myscalar;
            E_0(indkx,indkz,indwt) = E_0(indkx,indkz,indwt)*myscalar;
%             E_0(indkx,indkz,indwt) = E_0(indkx,indkz,indwt);
            
            X00 = Pbar0;
            % effect of base flow uncertainty at the level of alpha^2            
            X11 = lyap(A0+ii*wt*I2, A0', (1+ mum1)*Am1*X00 + (1+ mup1)*X00*Ap1');
            
            X02 = lyap(A0, (1 + mum1)*(Am1*X11' + X11*Am1') + (1 + mup1)*(Ap1*X11 + X11'*Ap1') ...
                         + sigma2p2*Ap1*X00*Ap1' + sigma2m2*Am1*X00*Am1');
           
            
            E_2(indkx,indkz,indwt) = real(trace(X02*Q));  % controlled flow with uncertain parameters
            
            
            % Total energy
            
            E(indkx, indkz,indwt) = E_0(indkx,indkz,indwt) + (alpha^2)*E_2(indkx,indkz,indwt);
            pre_E(indkx,indkz, indwt) = kx*kz*E(indkx,indkz,indwt);
            % The E2 correction without uncertainty
            pre_E2(indkx,indkz,indwt) = kx*kz*E_2(indkx,indkz,indwt);
            
            %% ===============================================================
            %% nuT_2 introduction
            %% obtain the second-order corrections to k, epsilon, and \nuT
            %% ===============================================================
            % velocity corelations in the flow with no control
            % the C operator mapping state to velocity field
            
            Cu = [ii*kx*D1, -ii*kz*I]/k2;
            Cv = [I, Z];
            Cw = [ii*kz*D1, ii*kx*I]/k2;
            
            % the D operator: D = \partail_y * C
            
            Du = [ii*kx*D2, -ii*kz*D1]/k2;
            Dv = [D1, Z];
            Dw = [ii*kz*D2, ii*kx*D1]/k2;
            
            uu_0(:,indkx,indkz,indwt) = real(diag(Cu*Pbar0*Cu'));
            uv_0(:,indkx,indkz,indwt) = real(diag(Cu*Pbar0*Cv'));
            uw_0(:,indkx,indkz,indwt) = real(diag(Cu*Pbar0*Cw'));
            vv_0(:,indkx,indkz,indwt) = real(diag(Cv*Pbar0*Cv'));
            vw_0(:,indkx,indkz,indwt) = real(diag(Cv*Pbar0*Cw'));
            ww_0(:,indkx,indkz,indwt) = real(diag(Cw*Pbar0*Cw'));

            uyvx_0(:,indkx,indkz,indwt) = real(-ii*kx*diag(Du*Pbar0*Cv'));
            uyuy_0(:,indkx,indkz,indwt) = real(diag(Du*Pbar0*Du'-R/2*(invIW/diag(1+nuT))*diag(E_turb_profile(:,indkx,indkz))*myscalar));
            vyvy_0(:,indkx,indkz,indwt) = real(diag(Dv*Pbar0*Dv'));
            vzwy_0(:,indkx,indkz,indwt) = real(ii*kz*diag(Cv*Pbar0*Dw'));
            wywy_0(:,indkx,indkz,indwt) = real(diag(Dw*Pbar0*Dw'-R/2*(invIW/diag(1+nuT))*diag(E_turb_profile(:,indkx,indkz))*myscalar));
            
            uxux_0(:,indkx,indkz,indwt) = kx^2*uu_0(:,indkx,indkz,indwt);
            wzwz_0(:,indkx,indkz,indwt) = kz^2*ww_0(:,indkx,indkz,indwt);
            uzwx_0(:,indkx,indkz,indwt) = kx*kz*uw_0(:,indkx,indkz,indwt);
            vxvx_0(:,indkx,indkz,indwt) = kx^2*vv_0(:,indkx,indkz,indwt);
            wxwx_0(:,indkx,indkz,indwt) = kx^2*ww_0(:,indkx,indkz,indwt);
            uzuz_0(:,indkx,indkz,indwt) = kz^2*uu_0(:,indkx,indkz,indwt);
            vzvz_0(:,indkx,indkz,indwt) = kz^2*vv_0(:,indkx,indkz,indwt);
            
            % 2nd-order correction to velocity corelations in the flow with control
            
            uu_2(:,indkx,indkz,indwt) = real(diag(Cu*X02*Cu'));
            uv_2(:,indkx,indkz,indwt) = real(diag(Cu*X02*Cv'));
            uw_2(:,indkx,indkz,indwt) = real(diag(Cu*X02*Cw'));
            vv_2(:,indkx,indkz,indwt) = real(diag(Cv*X02*Cv'));
            vw_2(:,indkx,indkz,indwt) = real(diag(Cv*X02*Cw'));
            ww_2(:,indkx,indkz,indwt) = real(diag(Cw*X02*Cw'));

            uyvx_2(:,indkx,indkz,indwt) = real(-ii*kx*diag(Du*X02*Cv'));
            uyuy_2(:,indkx,indkz,indwt) = real(diag(Du*X02*Du'));
            vyvy_2(:,indkx,indkz,indwt) = real(diag(Dv*X02*Dv'));
            vzwy_2(:,indkx,indkz,indwt) = real(ii*kz*diag(Cv*X02*Dw'));
            wywy_2(:,indkx,indkz,indwt) = real(diag(Dw*X02*Dw'));
            
            uxux_2(:,indkx,indkz,indwt) = kx^2*uu_2(:,indkx,indkz,indwt);
            wzwz_2(:,indkx,indkz,indwt) = kz^2*ww_2(:,indkx,indkz,indwt);
            uzwx_2(:,indkx,indkz,indwt) = kx*kz*uw_2(:,indkx,indkz,indwt);
            vxvx_2(:,indkx,indkz,indwt) = kx^2*vv_2(:,indkx,indkz,indwt);
            wxwx_2(:,indkx,indkz,indwt) = kx^2*ww_2(:,indkx,indkz,indwt);
            uzuz_2(:,indkx,indkz,indwt) = kz^2*uu_2(:,indkx,indkz,indwt);
            vzvz_2(:,indkx,indkz,indwt) = kz^2*vv_2(:,indkx,indkz,indwt);
            
            % k = k_0 + alpha^2 * k_2
            k_2(:,indkx,indkz,indwt) = (uu_2(:,indkx,indkz,indwt)+vv_2(:,indkx,indkz,indwt)+ww_2(:,indkx,indkz,indwt))/2;
            
            % epsilon = epsilon_0 + alpha^2 * epsilon_2
            epsilon_2(:,indkx,indkz,indwt) = (...
                2*(uxux_2(:,indkx,indkz,indwt) + vyvy_2(:,indkx,indkz,indwt) + wzwz_2(:,indkx,indkz,indwt) + ...
                uyvx_2(:,indkx,indkz,indwt) + uzwx_2(:,indkx,indkz,indwt) + vzwy_2(:,indkx,indkz,indwt)) + ...
                uyuy_2(:,indkx,indkz,indwt) + wywy_2(:,indkx,indkz,indwt) + vxvx_2(:,indkx,indkz,indwt) + ...
                wxwx_2(:,indkx,indkz,indwt) + uzuz_2(:,indkx,indkz,indwt) + vzvz_2(:,indkx,indkz,indwt) ...
                );
            
            % nu_T = nuT_0 + alpha^2 * nuT_2
            nuT_2(:,indkx,indkz,indwt) = diag(nuT)*(2*(diag(k_0RT)\k_2(:,indkx,indkz,indwt)) - diag(epsilon_0RT)\epsilon_2(:,indkx,indkz,indwt));
            
            % the 2nd-order correction to U if the pressure gradient is
            % kept constant and the bulk flux is allowed to change
            U2_nuT2(:,indkx,indkz,indwt) = LHS\(-(diag(nuT_2(:,indkx,indkz,indwt))*D2+diag(D1*nuT_2(:,indkx,indkz,indwt))*D1)*U0);
            
            % the 2nd-order correction to the bulk flux if the pressure 
            % gradient is kept constant and the bulk flux is allowed to change
            UB_nuT2(indkx,indkz,indwt) = diag(IW)'*U2_nuT2(:,indkx,indkz,indwt)/2;

            % the 2nd-order correction to pressure gradient for constant bulk flux
            % in our JFM 2011 paper, this is denoted by -\tau_{w2}
            % this is equal to the saved power
            % Psave = alpha^2 Px_nuT2
            Px_nuT2(indkx,indkz,indwt) = UB_nuT2(indkx,indkz,indwt)/(UB);
            
            % the 2nd-order correction to U if the pressure gradient is
            % adjusted to obtain constant bulk flux
            U2_UBconst_nuT2(:,indkx,indkz,indwt) = U2_nuT2(:,indkx,indkz,indwt) - U0*Px_nuT2(indkx,indkz,indwt);
            
            
            
            [indkx indkz indwt]
        end
    end
end

% integration scale for integrating the functions in log scale in kx and kz
logkxkz_scale = log(kxval(2)/kxval(1))*log(kzval(2)/kzval(1));

% premultiply by kx*kz for integration in log scale
for indkx = 1:kxgrd
    
    kx = kxval(indkx);
    
    for indkz = 1:kzgrd
        kz = kzval(indkz);
        k_2_pre(:,indkx,indkz,:)       = kx*kz*k_2(:,indkx,indkz,:);
        epsilon_2_pre(:,indkx,indkz,:) = kx*kz*epsilon_2(:,indkx,indkz,:);
        nuT_2_pre(:,indkx,indkz,:)     = kx*kz*nuT_2(:,indkx,indkz,:);
        U2_nuT2_pre(:,indkx,indkz,:)         = kx*kz*U2_nuT2(:,indkx,indkz,:);
        UB_nuT2_pre(indkx,indkz,:) = kx*kz*UB_nuT2(indkx,indkz,:);
        Px_nuT2_pre(indkx,indkz,:) = kx*kz*Px_nuT2(indkx,indkz,:);
        pre_nuT_2_avg(indkx,indkz) = kx*kz*sum(IW*nuT_2(:,indkx,indkz));
        pre_U2_nuT2_pre(indkx,indkz) = kx*kz*sum(IW*U2_nuT2(:,indkx,indkz));
    end
    
end
          
% integrate over horizontal wavenumbers
for indwt = 1:wtgrd
    
    k_2_avg(:,indwt)       = sum(sum(k_2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    epsilon_2_avg(:,indwt) = sum(sum(epsilon_2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    nuT_2_avg(:,indwt)     = sum(sum(nuT_2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    U2_nuT2_avg(:,indwt)          = sum(sum(U2_nuT2_pre(:,:,:,indwt),3),2)*logkxkz_scale;
    UB_nuT2_avg(indwt) = sum(sum(UB_nuT2_pre(:,:,indwt),2),1)*logkxkz_scale;
    Px_nuT2_avg(indwt) = sum(sum(Px_nuT2_pre(:,:,indwt),2),1)*logkxkz_scale;
    
end

rrr

            
     