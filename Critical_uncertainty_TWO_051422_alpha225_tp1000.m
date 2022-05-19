%% 05/14/2022
%% Written by Arya Abootorabi
%% Channel flow subject to randwt lower wall oscillations
% Optimal value for sigma2 

clear all
format long
format compact
clc

%%============================
% data for norm cwtputation %%
%%============================
% number of collocation points in y
N = 11;

% Reynolds numbers
R = 186;

% wall oscillaiton frequency
Tplusval = 1000;
wt = 2*pi*R/Tplusval;

% Oscillation amp
alpha = 2.25;

% imaginary unit
ii = sqrt(-1);

% identity and zero matrices
I = eye(N);
Z = zeros(N,N);
I2 = eye(2*N); 
ZZ = zeros(2*N, 2*N);

% parameters for uncertainty gamma_alpha, gamma_theta 

mu_alpha = 0;
mu_theta = 0;

% grid for kz - linear scale
% grid for kx, kz 
%kxlog = logspace(-4,0.48,51);
%kxval = [0 kxlog]; % important to include kx = 0!!!
kxval = 2.5;
%kzval = linspace(0.01,5,50);
kzval = 7;

% number of grid points in kx and kz plane
kxgrd = length(kxval);
kzgrd = length(kzval);

% =========================================================================
% Differentiation and integration matrices %
% =========================================================================
% Differentiation matrices
[yT,DM] = chebdif(N+2,3);
yvec=yT(2:end-1);

D1T = DM(:,:,1);
D2T = DM(:,:,2);
D3T = DM(:,:,3);

% implement hwtogeneous boundary coFnditions
D1 = DM(2:N+1,2:N+1,1);
D2 = DM(2:N+1,2:N+1,2);
D3 = DM(2:N+1,2:N+1,3);

% fourth derivative with clamped conditions
[y,D4]=cheb4c(N+2);

% integration weights
IWT=INTweights(N+2,2);
IW=IWT(2:end-1,2:end-1);
IW2 = kron(eye(2),IW);
IW3 = kron(eye(3),IW);
IWhalf = sqrtm(IW);
invIWhalf = IWhalf\I;
IW3half = sqrtm(IW3);
invIW3half = IW3half\IW3;

Id = eye(2);
Z2 = zeros(2,2);

% streamwise mean velocity
U = diag(1-yvec.^2);
Uy = diag(-2*yvec);
Uyy = diag(-2*ones(size(yvec)));

% Ws and Wc
fs = (1/2)*(1-yvec);
Wsbar = -inv(D2+((R*wt)^2)*inv(D2))*((R*wt)^2)*inv(D2)*fs;
Wcbar = inv(D2)*(R*wt)*(Wsbar+fs);
Ws = Wsbar + fs;
Wc = Wcbar;

% Derivatives of Ws and Wc
D1Ws = (D1*Wsbar)-((1/2)*ones(size(yvec)));  
D1Wc = D1*Wcbar;

D2Ws = D2*Wsbar;
D2Wc = D2*Wcbar;

Ws = diag(Ws);
Wc = diag(Wc);
D1Ws = diag(D1Ws);
D1Wc = diag(D1Wc);
D2Ws = diag(D2Ws);
D2Wc = diag(D2Wc); 


% When W is written with exponential terms
Wp1 = Wc - (ii*Ws);
Wm1 = Wc + (ii*Ws);

% second derivatives of Wp1 and Wm1
D2Wp1 = D2Wc - (ii*D2Ws);
D2Wm1 = D2Wc + (ii*D2Ws);

D1Wp1 = D1Wc - (ii*D1Ws);
D1Wm1 = D1Wc + (ii*D1Ws);

%%=========================================================================
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
% UB = sum(IW'*U0)/2;

% Wp
Wp1bar = (LHS/R-ii*wt*I)\(wt*ones(N,1));
Wp1ybar = D1*Wp1bar;
Wp1yybar = D2*Wp1bar;
Wp1 = Wp1bar + (-ii)*ones(N,1);
Wp1y = Wp1ybar;
Wp1yy = Wp1yybar;

% Wm
Wm1bar = (LHS/R+ii*wt*I)\(wt*ones(N,1));
Wm1ybar = D1*Wm1bar;
Wm1yybar = D2*Wm1bar;
Wm1 = Wm1bar + (ii)*ones(N,1);
Wm1y = Wm1ybar;
Wm1yy = Wm1yybar;

% number of harmonics
m = 5;

% matrix Emat
Emat = -ii* m* wt* I2;

for m = -(m-1):m
    Emat = blkdiag(Emat,ii*m*wt*I2);
end

% grid for sigma2_alpha
sigma2alphaval = linspace(1e-6,5e-2,10);
% sigma2alphaval = 1e-2;

% power iteration parameters
lamgrd = 30;
lam_residual = 1e-4;
sigma2alphagrd = length(sigma2alphaval);

% Initializing the vector for sigma2theta
sigma2thetaval = zeros(1,sigma2alphagrd);

tic

for indkx = 1:kxgrd
    
    kx = kxval(indkx);
    
    for indkz = 1:kzgrd
        
        kz = kzval(indkz);
        
        %%=====================
        %% k2 := kx^2 + kz^2 %%
        %%=====================
        k2 =kx^2 + kz^2;
        k4 = k2*k2;
        
        %%=============
        %% Laplacian %%
        %%=============
        Delta = D2 - k2*I;
        Delta2 = D4 - 2*k2*D2 + k4*I;
        
         % Q matrix
        QvT = (1/k2)*D1T'*IWT*D1T + IWT;
        QetaT = (1/k2)*IWT;

        Qv = QvT(2:N+1,2:N+1);
        Qeta = QetaT(2:N+1,2:N+1);
        Q = [Qv, Z; Z, Qeta];
        
        Im = eye(2*m+1,2*m+1);
        bigQ = kron(Im, Q);
        
        % square roots of Qv and Qeta
        Qv_half = sqrtm(Qv);
        Qeta_half = sqrtm(Qeta);

        Qhalf = [Qv_half, Z; Z, Qeta_half];
        bigQhalf = kron(Im, Qhalf);
        
        invQhalf = [Qv_half\I, Z; Z, Qeta_half\I];
        biginvQhalf = kron(Im, invQhalf);
  
        % inverse of Q
        invQv = Qv\I;
        invQeta = Qeta\I;
        invQ = [invQv, Z; Z, invQeta];
        
        biginvQ = kron(Im, invQ);
        
        %% Matrix A

        % Matrix A0 for channel flow with no wall oscillations
        A11 = Delta\(1/R*(diag(1+nuT))*Delta2 + 1/R*diag(nuTyy)*(D2+k2*I) + 2/R*diag(nuTy)*(D3 - k2*D1) + ...
            ii*kx*(diag(U0yy) - diag(U0)*Delta));
        A12 = Z;
        A21 = -ii*kz*diag(U0y);
        A22 = 1/R*diag(1+nuT)*Delta + 1/R*diag(nuTy)*D1 - ii*kx*diag(U0);
        
        A0 = [A11, A12; A21, A22];

        
        %% Ap1 and Am1
        
        % Am1 operator (order of alpha^1)
        
        A11 = Delta\(ii*kz*(diag(Wm1yy) - diag(Wm1)*Delta));
        A12 = Z;
        A21 = ii*kx*diag(Wm1y);
        A22 = -ii*kz*diag(Wm1);
        
        Am1 = [A11, A12; A21, A22];

        % Ap1 operator (order of alpha^1)
        
        A11 = Delta\(ii*kz*(diag(Wp1yy) - diag(Wp1)*Delta));
        A12 = Z;
        A21 = ii*kx*diag(Wp1y);
        A22 = -ii*kz*diag(Wp1);
        
        Ap1 = [A11, A12; A21, A22];
        
        %%  The crticial uncertainty
        
        for indsigmaalpha = 1:sigma2alphagrd
            
            sigma2alpha = sigma2alphaval(indsigmaalpha);
            
            % bounds for sigma2
            sigma2theta_lower = 0;
            sigma2theta_upper = 1;
            interlen = sigma2theta_upper - sigma2theta_lower;
            
            while interlen > lam_residual 
                sigma2theta_mid = (sigma2theta_lower + sigma2theta_upper)/2;
                mup1 = (1 + mu_alpha)*exp(ii*mu_theta + sigma2theta_mid/2) - 1;
                mum1 = (1 + mu_alpha)*exp(-ii*mu_theta + sigma2theta_mid/2) - 1;

                sigma2p1 = (exp(2*ii*mu_theta + sigma2theta_mid))*(sigma2alpha*exp(sigma2theta_mid) + ((1 + mu_alpha)^2)*(exp(sigma2theta_mid) - 1));
                sigma2m1 = (exp(-2*ii*mu_theta + sigma2theta_mid))*(sigma2alpha*exp(sigma2theta_mid) + ((1 + mu_alpha)^2)*(exp(sigma2theta_mid) - 1));

                % truncated toeplitz A matrix; nwtinal
                Abar = full(blktridiag(A0,alpha*(1+mup1)*Ap1,alpha*(1+mum1)*Am1,2*m+1))- Emat;

                eigAbar = real(eig(Abar));

                    if(~sum( eigAbar > 0 ))

                        % B0 matrix ;( [Am1, Ap1] )
                        bigAp1 = full(blktridiag(ZZ,alpha*Ap1,ZZ,2*m+1));
                        bigAm1 = full(blktridiag(ZZ,ZZ,alpha*Am1,2*m+1));

                        B0 = [bigAm1, bigAp1];

                        % initialization of randwt hermitian matrix P
                        P0tmp = rand(2*N*(2*m+1), 2*N*(2*m+1));
                        P0tmp = (P0tmp + P0tmp')/2;
                        eigP0min = min(real(eig(P0tmp)));
                        if eigP0min < 0
                            P0tmp = 2*abs(eigP0min)*eye(2*N*(2*m+1)) + P0tmp;
                        end
                        P0tmp = P0tmp/ norm(P0tmp,'fro');

                        P0 = [P0tmp P0tmp; P0tmp P0tmp];            
                        Xker = lyap(Abar, B0*P0*B0');
                        Xker = (Xker + Xker')/2;

                        P1tmp = [sigma2p1*Xker, zeros(2*N*(2*m+1)); zeros(2*N*(2*m+1)), sigma2p1*Xker];

                        % Power iteration for spectral radius of \cL (loop-gain) operator
                        lam = [];
                        
                            for indlam = 1:lamgrd
                                P0 = P1tmp/norm(P1tmp,'fro');
                                Xker = lyap(Abar, B0*P0*B0');
                                Xker = (Xker + Xker')/2;
                                P1tmp = [sigma2m1*Xker, zeros(2*N*(2*m+1)); zeros(2*N*(2*m+1)), sigma2p1*Xker];
                                lam(indlam) = trace(P0'*P1tmp);
                                r = P1tmp - P0*lam(indlam);   % relative error
                                normr = norm(r,'fro')/norm(P1tmp,'fro');
                                if mod(indlam,10) == 0
                                   [lam']
                                   pause(0.1);
                                end
                                if normr < lam_residual
                                    break
                                end
                            end

                            %lam(end) = 1.00009999
                            diff_lam = lam(end) - 1;

                            if diff_lam < lam_residual % we can afford to be brave and increase sigma
                                sigma2theta_lower = sigma2theta_mid;
                            else % we need to be conservative and lower sigma
                                sigma2theta_upper = sigma2theta_mid;
                            end
                    else
                        sigma2theta_upper = sigma2theta_mid;
                    end
                interlen = sigma2theta_upper - sigma2theta_lower;    
                disp(['indsigmaalpha = ',num2str(indsigmaalpha),'; interlen1 = ',num2str(interlen)]);
            end
            sigma2thetaval(indsigmaalpha) = sigma2theta_mid;
        end            
            
         [indkx,indkz]
    end
end

rrr

str1 = 'sigmas';
str2 = num2str(indkx);
str3 = num2str(indkz);
name = strcat(str1,'_indkx',str2,'_indkz',str3);
save(name,'sigma2alphaval','sigma2thetaval')
time = toc;

