function [pN1_b,pN2_b,P]=SDFFTPred2(DataSheet,fNo,fNg,s)
% Analytically find the find steady-state distribution P(i,j,k). 
% For a fixed block k, P(i,j) gives the probability of observing i red and 
% j blue agents. 
% DataSheet: the initial steady-state data from Schelling1
% fNo: total number of red agent after the regional demographic change
% fNg: total number of blue agent after the regional demographic change
% s: available cells in each block, calculated in Schelling1 
% pN1_b(i): Mean number of red agents in i-th block 
% pN2_b(i): Mean number of blue agents in i-th block 
%% Find potential and frustration (DFFT)
    [h] = get_h(DataSheet,s);
    [~,~,~,f,V1,V2,~,~,~,~,~,~,~,~] = TC_DFFT_Y_2(h,8,1,s);
    f=RemovePlane(f);
    V1=V1-mean(V1);
    V2=V2-mean(V2);
    %% Inputs
    % change potential into the right dimension
    V1=reshape(V1,[1,1,length(V1)]);
    V2=reshape(V2,[1,1,length(V2)]);
    % Determine spots where there is no data and replace numbers in those spots
    % of f with NaN
    % h=hist3([occ1(:),occ2(:)],...
    %     'Ctrs',{0:max(occ1(:)) 0:max(occ2(:))},'CdataMode','auto');
    % f(h==0)=NaN;
    totN1=fNo;
    totN2=fNg;
    mu1_s=0; % Adjust Initial guess so that Newton's method converges 
    mu2_s=0; % Adjust Initial guess so that Newton's method converges 
    eps=1e-5;
    steps=10;

    %% Calculate Mu's
    [Mu,f_Mu]=FindMu2(totN1,totN2,V1,V2,f,mu1_s,mu2_s,eps,steps,s)

    %% Find expected number of people for each location
    % Define useful matrices:
    ROI=size(f);
    N1_mat=repmat((0:ROI(1)-1)',[1,ROI(2)]);
    N2_mat=repmat((0:ROI(2)-1),[ROI(1),1]);
    %InvFactProd=1./(factorial(N1_mat).*factorial(N2_mat));
    Ne=s-N1_mat-N2_mat;
    Ne(Ne<0)=NaN;
    logInvFactProd=-gammaln(N1_mat+1)-gammaln(N2_mat+1)-gammaln(Ne+1); 
    
    % Find unnormalized probability distribution
    logP=logInvFactProd-(V1-Mu(1)).*N1_mat-(V2-Mu(2)).*N2_mat-f;
    P=exp(logP);
    % convert NaN to 0 (Since these values are not observed, they 
    % should have probability~0 ??? Not really. Effectively saying that
    % f are Inf at these spots, clearly not the case all the time.)
    P(isnan(P))=0;
    % Normalize P
    P=P./sum(sum(P,1),2);
    % Find expected number of N1 and N2 for this set of Mu:
    pN1_b=reshape(sum((0:ROI(1)-1)'.*sum(P,2),1),[],1);
    pN2_b=reshape(sum((0:ROI(2)-1).*sum(P,1),2),[],1);
end


%% Subfunctions 
function [h] = get_h(DataSheet,s)
%Get the joint histograms
    CellSize=size(DataSheet,1);
    h=zeros(s+1,s+1,CellSize^2); 
    k=0;
    [x,y]=meshgrid(0:s,0:s);
    tf=x+y>s;
    for i=1:CellSize %for each bin
        for j=1:CellSize
            k=k+1;
            Data=DataSheet{i,j};
            h(:,:,k)=hist3(Data,'Ctrs',{0:s 0:s},'CdataMode','auto');
        end
    end
end

function [Mu,f_Mu] = FindMu2(totN1,totN2,V1,V2,f,mu1_s,mu2_s,eps,steps,s)
% Find optimized chemical potentials mu1 and mu2 to match totN1 and totN2, 
% Mu=[mu1;mu2], f_Mu gives the deviation of current prediction to totN1/N2

% V1 and V2 are potentials of dimension 1*1*Nb, where Nb is the number of
% different locations

% F is frustration of dimension (maxN1+1)*(maxN2+1). Its entry is NaN if
% there is no data for that entry for its estimation

% mu1_s and mus2_s are starting points for the search

% eps is the displacment used in (mu1,mu2) space to approximate 
% its derivative 
    BAD=false;
    % Define useful matrices:
    ROI=size(f);
    N1_mat=repmat((0:ROI(1)-1)',[1,ROI(2)]);
    N2_mat=repmat((0:ROI(2)-1),[ROI(1),1]);

    Ne=s-N1_mat-N2_mat;
    Ne(Ne<0)=NaN;
    logInvFactProd=-gammaln(N1_mat+1)-gammaln(N2_mat+1)-gammaln(Ne+1); 
    % Start Newton's Method:
    Mu=[mu1_s;mu2_s];
    
    for k=1:steps
      
        % Evaluate function func at Mu. We're trying to find the root of
        % func using Newton's Method
        f_Mu=func2(Mu,V1,V2,f,totN1,totN2,logInvFactProd,N1_mat,N2_mat,ROI);
        % Evaluate function at Mu+[eps;0]
        Mu1=Mu+[eps;0];
        f_Mu1=func2(Mu1,V1,V2,f,totN1,totN2,logInvFactProd,N1_mat,N2_mat,ROI);
        % Evaluate function at Mu+[0;eps]
        Mu2=Mu+[0;eps];
        f_Mu2=func2(Mu2,V1,V2,f,totN1,totN2,logInvFactProd,N1_mat,N2_mat,ROI);      
        
        % Estimate Jacobian of f at Mu:
        Df=[(f_Mu1-f_Mu)/eps,(f_Mu2-f_Mu)/eps]';
        % Check if [Df] is invertible
        if det(Df)==0
            fprintf('[Df] not invertible at step %d\n.',k)
            BAD=true;
            break
        end
        % Find next Mu using Newton's method:
        Mu=Mu-Df\f_Mu;
    end
    if BAD
        Mu=NaN;
    end
end


function [f_Mu]=func2(Mu,V1,V2,f,totN1,totN2,logInvFactProd,N1_mat,N2_mat,ROI)
% Returns [predN1-totN1;predN2-totN2] for the set of mu's given in Mu
    % Find unnormalized probability distribution
    logP=logInvFactProd-(V1-Mu(1)).*N1_mat-(V2-Mu(2)).*N2_mat-f;
    P=exp(logP);
    % convert NaN to 0 (Since these values are not observed, they 
    % should have probability~0 ??? Not really. Effectively saying that
    % f are Inf at these spots, clearly not the case all the time.)
    P(isnan(P))=0;
    % Normalize P
    P=P./sum(sum(P,1),2);
    % Find expected number of N1 and N2 for this set of Mu:
    predN1=sum(sum((0:ROI(1)-1)'.*sum(P,2),1),3);
    predN2=sum(sum((0:ROI(2)-1).*sum(P,1),2),3);
    % Output of the function 
    f_Mu=[predN1-totN1;predN2-totN2];
end

function [srho]=RemovePlane(rho)
% Remove plane from a surface 
    d=size(rho);
    x=[];
    y=[];
    z=[];
    k=0;
    for i=1:d(1)
        for j=1:d(2)
            if ~isnan(rho(i,j))
                x=[x;i];
                y=[y;j];
                z=[z;rho(i,j)];
            end
        end
    end 
    f=fit([x,y],z,'poly11');

    % Convert f to matrix
    Fit=zeros(d);
    for i=1:d(1)
        for j=1:d(2)
            Fit(i,j)=f(i,j);
        end
    end   
    srho=rho-Fit;
end



