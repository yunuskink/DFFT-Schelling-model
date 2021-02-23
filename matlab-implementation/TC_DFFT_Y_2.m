function [P,rho,sigma_rho,f,V1,V2,P_pred,chi_sq,lglik,f_params,f_poly,mx_bin1,mx_bin2,set_gauge] = TC_DFFT_Y_2(h,order_f,tau,s)

%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%
%h-> Two component histogram, can be calculated using "get_hist" from
%density observations.
%order_f -> Highest order polynomial to use to fit to f. A value of 8 is as
%high as you ever need to go. If there is a warning about a gauge error,
%decreasing this and/or removing outliers and reducing the limits of h will
%get rid of that warning, but is not essential to do so.
%tau -> Autocorrelation constant for dealing with highly correlated data.
%Used for calculating uncertainties. Unnecessary, default value is 1.

%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%
%P -> Probability distribution from the data calculated from histogram, no
%modelling used.
%rho-> -log(N1!*N2!*P). Parameters are fit to this variable.
%sigma_rho -> uncertainty in above variable.
%f -> Frustration for each observed density value. Size
%[N1_max+1,N2_max+1]. NaN value are for values where h==0.
%V1 -> Vexation for component one for each bin.
%V2 -> Vexation for component two for each bin.
%P_pred -> Probability distribution from f, V1, and V2.
%chi_sq -> Result of goodness of fit test.
%lglik -> Log-likelihood value
%Rest of params are needed to fit into MLE_poly to get the MLE estimate of
%the parameters and eventually the uncertainty.


%%%%%%%%%%%%%%DATASET PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_bins = size(h,3);
N1_max = size(h,1)-1; 
N2_max = size(h,2)-1;
if nargin<3
    tau = 1;
end
set_gauge = 1;
gauge_wght = 20;

%%%%%%%%%%%%%%%%%%%%%P = P(N_1,N_2,b)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = zeros(size(h));
N_obs = zeros(1,size(h,3));
for i=1:size(P,3)
    N_obs(i) = (sum(sum(h(:,:,i),1),2));
    P(:,:,i) = h(:,:,i)./N_obs(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%SIGMA P%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_p = zeros(size(P));
sigma_rho = zeros(size(P));
for i=1:size(P,3)  %for each bin
    observations = h(:,:,i)/tau; %Independent obervations of each number of flies in bin i.
    sigma_p(:,:,i) = sqrt(observations.*(N_obs(i)-observations)/(N_obs(i)^2*(N_obs(i)+1))); %Equation (43)
    sigma_rho(:,:,i) = sqrt(psi(1,(h(:,:,i))/tau) - psi(1,(N_obs(i)/tau)));
end

%%%%%%%%%%%%%%%%%%%%%%BUILDING SET OF ORTHOGONAL POLYNOMIALS%%%%%%%%%%%%%%%
N1_mat = repmat([0:N1_max]',[1 N2_max+1 N_bins]);
N2_mat = repmat([0:N2_max],[N1_max+1 1 N_bins]);
x1_mat = (N1_mat(:,:,1)-N1_max/2)/(N1_max/2);
x2_mat = (N2_mat(:,:,1)-N2_max/2)/(N2_max/2);
N_poly = sum([2:order_f]+1); %INCLUDES 0th ORDER
f_poly = zeros(size(N1_mat,1),size(N1_mat,2),N_poly);
TC_orders = zeros(N_poly,2);

poly_ind = 1;
for i=2:order_f
    for j=0:i %NO 0th ORDER %for j=0:i %INCLUDES 0th ORDER
        TC_orders(poly_ind,1) = j;
        TC_orders(poly_ind,2) = i-TC_orders(poly_ind,1);
        f_poly(:,:,poly_ind) = polyval(LegendrePoly(TC_orders(poly_ind,1)), x1_mat).*polyval(LegendrePoly(TC_orders(poly_ind,2)), x2_mat);
        poly_ind = poly_ind+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATING RHO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaln_mat = gammaln(N1_mat+1) + gammaln(N2_mat+1);
[x,y]=meshgrid(0:s);
Ne=s-x-y;
Ne(Ne<0)=NaN;
Ne=Ne(1:(N1_max+1),1:(N2_max+1));
LogP=log(P);
LogP(LogP==-Inf)=NaN;
rho = -LogP - gammaln_mat -gammaln(Ne+1);

%%%%%%%%%%%%%%%%%BUILDING LEAST SQUARES MATRICES%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve LHS = RHS*parameters. LHS contains the data, and RHS contains the
%coefficients of the parameters of d/d theta for each parameter, theta.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LHS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS = zeros(N_bins*3+N_poly,1);
%%%%%%%%%%%%%%%%z_b%%%%%%%%%%%%%%%%%
for i=1:N_bins
    LHS(i) = nansum(nansum(rho(:,:,i)./(sigma_rho(:,:,i).^2)));
end
%%%%%%%%%%%%%%%%v_b1%%%%%%%%%%%%%%%%%
for i=1:N_bins
    LHS(N_bins+i) = nansum(nansum(N1_mat(:,:,i).*rho(:,:,i)./(sigma_rho(:,:,i).^2)));
end
%%%%%%%%%%%%%%%%v_b2%%%%%%%%%%%%%%%%%
for i=1:N_bins
    LHS(2*N_bins+i) = nansum(nansum(N2_mat(:,:,i).*rho(:,:,i)./(sigma_rho(:,:,i).^2)));
end

%%%%%%%%%%%%%%%%alpha_i%%%%%%%%%%%%%%%%%
for j=1:N_poly
    LHS(3*N_bins+j) = nansum(nansum(nansum(repmat(f_poly(:,:,j),[1 1 N_bins]).*rho./(sigma_rho.^2))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RHS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RHS = zeros(N_bins*3+N_poly);
%%%%%%%%%%%%%%%%z_b%%%%%%%%%%%%%%%%%
for i=1:N_bins %For each derivative
    for k=1:N_bins %z_b
        if i==k
            RHS(i,k) = sum(sum(1./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_bins %v_b1
        if i==k
            RHS(i,k+N_bins) = sum(sum(N1_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_bins %v_b2
        if i==k
            RHS(i,k+2*N_bins) = sum(sum(N2_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_poly %alpha_i
        RHS(i,k+3*N_bins) = sum(sum(f_poly(:,:,k)./(sigma_rho(:,:,i).^2)));
    end
end

%%%%%%%%%%%%%%%%v_b1%%%%%%%%%%%%%%%%%
for i=1:N_bins
    for k=1:N_bins %z_b
        if i==k
            RHS(N_bins+i,k) = sum(sum(N1_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_bins %v_b1
        if i==k
            RHS(N_bins+i,k+N_bins) = sum(sum(N1_mat(:,:,i).*N1_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_bins %v_b2
        if i==k
            RHS(N_bins+i,k+2*N_bins) = sum(sum(N1_mat(:,:,i).*N2_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
        end
    end
    for k=1:N_poly %alpha_i
        RHS(N_bins+i,k+3*N_bins) = sum(sum(N1_mat(:,:,i).*f_poly(:,:,k)./(sigma_rho(:,:,i).^2)));
    end
end

%%%%%%%%%%%%%%%%v_b2%%%%%%%%%%%%%%%%%
for i=1:N_bins
    for k=1:N_bins %z_b
        if i==k
            RHS(2*N_bins+i,k) = sum(sum(N2_mat(:,:,i)./(sigma_rho(:,:,i).^2)));
        end
    end
    for k=1:N_bins %v_b1
        if i==k
            RHS(2*N_bins+i,k+N_bins) = sum(sum(N2_mat(:,:,i).*N1_mat(:,:,i)./(sigma_rho(:,:,i).^2)));
        end
    end
    for k=1:N_bins %v_b2
        if i==k
            RHS(2*N_bins+i,k+2*N_bins) = sum(sum(N2_mat(:,:,i).*N2_mat(:,:,i)./(sigma_rho(:,:,i).^2)));
        end
    end
    for k=1:N_poly %alpha_i
        RHS(2*N_bins+i,k+3*N_bins) = sum(sum(N2_mat(:,:,i).*f_poly(:,:,k)./(sigma_rho(:,:,i).^2)));
    end
end

%%%%%%%%%%%%%%%%alpha_i%%%%%%%%%%%%%%%%%
for i=1:N_poly
    for k=1:N_bins %z_b
        RHS(3*N_bins+i,k) = sum(sum(f_poly(:,:,i)./(sigma_rho(:,:,k).^2)));
    end
    for k=1:N_bins %v_b1        
        RHS(3*N_bins+i,k+N_bins) = sum(sum(f_poly(:,:,i).*N1_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
    end
    for k=1:N_bins %v_b2
        RHS(3*N_bins+i,k+2*N_bins) = sum(sum(f_poly(:,:,i).*N2_mat(:,:,k)./(sigma_rho(:,:,k).^2)));
    end
    for k=1:N_poly %alpha_i
        for m=1:N_bins
            RHS(3*N_bins+i,k+3*N_bins) = RHS(3*N_bins+i,k+3*N_bins) + sum(sum(f_poly(:,:,i).*f_poly(:,:,k)./(sigma_rho(:,:,m).^2)));
        end
    end
end

%%%%%%%%%%%%%%%%GAUGE FIXING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last equations to fix the gauge.
%mx_bin1 = find(sum(occ1,1)==max(sum(occ1,1)));mx_bin1 = mx_bin1(1);
%mx_bin2 = find(sum(occ2,1)==max(sum(occ2,1)));mx_bin2 = mx_bin2(1);
mx_bin1 = 1;
mx_bin2 = 1;

if set_gauge    
    %Slope gauges
    LHS = [LHS;0;0];
    vect1 = zeros(1,size(RHS,2));
    vect1(1*N_bins+mx_bin1) = gauge_wght;
    vect2 = zeros(1,size(RHS,2));
    vect2(2*N_bins+mx_bin2) = gauge_wght;
    RHS = [RHS;vect1;vect2];
    %Constant shift gauge
    %LHS = [LHS;N_poly];
    %vect3 = zeros(1,size(RHS,2));
    %vect3(3*N_bins+1:end) = squeeze(f_poly(N1_max+1,N2_max+1,:))*gauge_wght;
    %RHS = [RHS;vect3];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOLVING LINEAR PROBLEM%%%%%%%%%%%%%%%%%%%%%
sol = RHS\LHS;
z_b = sol(1:N_bins);
V1 = sol(N_bins+1:2*N_bins);
V2 = sol(2*N_bins+1:3*N_bins);
f_params = sol(3*N_bins+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%BUILDING f FROM POLYNOMIALS%%%%%%%%%%%%%%%%%%%%
f = zeros(size(f_poly,1),size(f_poly,2));
observed_densities = sum(h,3)>=1;
for i=1:size(f_poly,3)
    f = f+f_poly(:,:,i).*f_params(i);
end
f(~observed_densities) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%CHECKING NORMALIZATION OF P%%%%%%%%%%%%%%%%%%%
P_sol = zeros(size(rho));
for i=1:size(P_sol,3)
    P_sol(:,:,i) = exp(-(z_b(i)+N1_mat(:,:,i).*V1(i)+N2_mat(:,:,i).*V2(i)+f+gammaln_mat(:,:,i)));
    P_sol(:,:,i) = (h(:,:,i)>=1).*P_sol(:,:,i);
end
histogram(nansum(nansum(P_sol,1),2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATING PROPERLY NORMALIZED P_pred%%%%%%%%%%%
P_pred = zeros(size(rho));
for i=1:size(P_pred,3)
    P_pred(:,:,i) = exp(-(N1_mat(:,:,i).*V1(i)+N2_mat(:,:,i).*V2(i)+f+gammaln_mat(:,:,i)));
    P_pred(:,:,i) = P_pred(:,:,i)./nansum(nansum(P_pred(:,:,i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%GOODNESS OF FIT TESTS%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_model = zeros(size(rho));
for i=1:size(rho_model,3)
    rho_model(:,:,i) = z_b(i)+N1_mat(:,:,i).*V1(i)+N2_mat(:,:,i).*V2(i)+f;
end
chi_sq_mat = sqrt((rho-rho_model).^2)./sigma_rho;
chi_sq = nansum(chi_sq_mat(:));

lglik = nansum(nansum(nansum(h.*log(P_pred))));


% LegendrePoly.m by David Terr, Raytheon, 5-10-04

% Given nonnegative integer n, compute the 
% Legendre polynomial P_n. Return the result as a vector whose mth
% element is the coefficient of x^(n+1-m).
% polyval(LegendrePoly(n),x) evaluates P_n(x).


function pk = LegendrePoly(n)

if n==0 
    pk = 1;
elseif n==1
    pk = [1 0]';
else
    
    pkm2 = zeros(n+1,1);
    pkm2(n+1) = 1;
    pkm1 = zeros(n+1,1);
    pkm1(n) = 1;

    for k=2:n
        
        pk = zeros(n+1,1);

        for e=n-k+1:2:n
            pk(e) = (2*k-1)*pkm1(e+1) + (1-k)*pkm2(e);
        end
        
        pk(n+1) = pk(n+1) + (1-k)*pkm2(n+1);
        pk = pk/k;
        
        if k<n
            pkm2 = pkm1;
            pkm1 = pk;
        end
        
    end
    
end