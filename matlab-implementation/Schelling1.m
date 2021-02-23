% Obtain Steady-state Ensemble; Create Altered Emsemble


% OUTPUTS:
% DataSheet: steady-state ensemble. DataSheet{i,j}(k,1) gives number of red
% in {i,j}-th block for k-th ensemble.DataSheet{i,j}(k,2) gives number of
% blue in {i,j}-th block for k-th city in the ensemble.
% Plist, Elist: Altered ensemble, to be used in Schelling2 and predictions.
% Plist(i,1/2/3,k) gives the x position/y position/tag of the i=th agent in
% k-th city in the ensemble. tag=0 means red; tag=1 means blue. The 
% additional information of location is only needed for Schelling2 but not 
% predictions.
% Omove: Observed number of moves in the city (to find time scale)
%% Parameters:
iNo=1000; % Initial number of Red 
iNg=1000; % Initial number of Blue
rr=350; % number of red replaced by blue
re=0; % number of empty replaced by blue


SampleSize=2500; % Size of emsemble - smaller for speed for demonstration purposes
%SampleSize=50000; % Size of emsemble used in ArXiv paper
SampleStep=1000; % Number of steps to sample over time to get steady-state ensemble
Order=60; % Schelling lattice grid size (Order*Order possible positions) 
CellSize=5; % DFFT coarse-grain block grid size: a number that divides Order
Nblock=CellSize^2; % Number of blocks in the city
s=(Order/CellSize)^2; % Number of cells in each block  
r=1; % neighbor radius: only people within r/Order are considered neighbors



%% Main 

% Initialize
OMove=0;
No=iNo;
Ng=iNg; 
N=No+Ng;
ep=1/(2*Order); % correct for float point error when comparing sizes
tag=[zeros(No,1);ones(Ng,1)]; % Tag for agent color. R:0, B:1
InitialP1=zeros(No,2,SampleSize);
InitialP2=zeros(Ng,2,SampleSize);
InitialE=zeros(Order^2-No-Ng,2,SampleSize);
DataSheet=cell(CellSize);
for m=1:CellSize
    for n=1:CellSize
        DataSheet{m,n}=zeros(SampleSize,2);
    end
end



% Randomly initialize positions: Each row contains to x,y position of one
% person/cell
[Pp,Empty]=GenerateRandom(Order,No,Ng);
P=[Pp,tag];
LEmpty=size(Empty,1);

% Plot(P,No,Order)
% Initial Run
for i=1:100000
    % Randomly pick a pair of indices 
    pi=ceil(rand*N); % random person
    qi=N+ceil(rand*LEmpty); % random empty spot
    % Find total utility/potential before switch
    U=FindUtility(pi,P,Order,r);
    Pi=Potential(pi,P);
    % Find P after switch
    LOC=[P(:,1:2);Empty];
    p=LOC(pi,:);
    q=LOC(qi,:);
    LOC(pi,:)=q;
    LOC(qi,:)=p;
    sP=[LOC(1:N,:),tag];
    sEmpty=LOC(N+1:end,:);
    V=FindUtility(pi,sP,Order,r);
    Pf=Potential(pi,sP);
    % Move with probability
    if rand<1/(1+exp((U+Pi)-(V+Pf)))
        P=sP;
        Empty=sEmpty;
    end
end

% Run
for i=1:SampleSize
    if mod(i,100)==0
        disp(i)
    end
    for k=1:SampleStep
        % Randomly pick a pair of indices 
        pi=ceil(rand*N); % random person
        qi=N+ceil(rand*LEmpty); % random empty spot
        % Find total utility/potential before switch
        U=FindUtility(pi,P,Order,r);
        Pi=Potential(pi,P);
        % Find P after switch
        LOC=[P(:,1:2);Empty];
        p=LOC(pi,:); % find and switch coordinates
        q=LOC(qi,:);
        LOC(pi,:)=q;
        LOC(qi,:)=p;
        sP=[LOC(1:N,:),tag];
        sEmpty=LOC(N+1:end,:);
        V=FindUtility(pi,sP,Order,r);
        Pf=Potential(pi,sP);
        % Move with probability
        if rand<1/(1+exp((U+Pi)-(V+Pf)))
            P=sP;
            Empty=sEmpty;
            OMove=OMove+1;
        end
    end
    % Sort result into DataSheet
    InitialP1(:,:,i)=P(1:No,1:2);
    InitialP2(:,:,i)=P(No+1:end,1:2);
    InitialE(:,:,i)=Empty;
    % Sort result into DataSheet
    for m=1:CellSize % Row
        for n=1:CellSize % Column
            % Find number of orange and green people:
            DataSheet{m,n}(i,1)=sum((P(:,3)==0).*(P(:,1)>=(n-1)/CellSize-ep)...
                .*(P(:,1)<(n/CellSize-ep)).*(P(:,2)>=(1-m/CellSize-ep))...
                .*(P(:,2)<(1-(m-1)/CellSize-ep)));
            DataSheet{m,n}(i,2)=sum((P(:,3)==1).*(P(:,1)>=(n-1)/CellSize-ep)...
                .*(P(:,1)<(n/CellSize-ep)).*(P(:,2)>=(1-m/CellSize-ep))...
                .*(P(:,2)<(1-(m-1)/CellSize-ep)));
        end
    end
end

%% Sample plot showing the steady-state Schelling city
figure
Plot(P,No,Order)

%% Alter ensemble
N=iNo+iNg+re;
[Plist,Elist]=GenerateInitial7(InitialP1,InitialP2,InitialE,rr,re,iNo,iNg,Order,r,N);

%% Plot altered ensemble
figure
Plot(Plist(:,:,end),iNo-rr,Order)

%% SUB Functions 
function [P,Empty] = GenerateRandom(Order,No,Ng)
% generate a random 
    [X,Y]=meshgrid(0:(Order-1));
    x=reshape(X,[],1);
    y=reshape(Y,[],1);
    L=length(x);
    i=datasample(1:L,No+Ng,'Replace',false);
    P=[x(i),y(i)]/Order;
    tfEmpty=true(L,1);
    tfEmpty(i)=false;
    Empty=[x(tfEmpty),y(tfEmpty)]/Order;
end

function [u]=FindUtility(pi,P,Order,r)
% Find utility of agent with index p
    p=P(pi,1:2);
    og=P(pi,3);
    [NO,NG]=FindNeighborN(p,P,Order,r,og); % Find Number of neighbors of p
    u=Utility(NO,NG,og);
end

function [NO,NG] = FindNeighborN(p,P,Order,r,og)
% Convert coordinates of P into indices of a matrix (flipped, extended), 
% then find neighbors around p for each color
    O=P(P(:,3)==0,1:2);
    O=round(O*Order+1+r);
    p=round(p*Order+1+r);
    Lindex=O(:,1)+(O(:,2)-1)*(Order+2*r);
    MO=zeros(Order+2*r);
    MO(Lindex)=1;
    G=P(P(:,3)==1,1:2);
    G=round(G*Order+1+r);
    Lindex=G(:,1)+(G(:,2)-1)*(Order+2*r);
    MO(Lindex)=1i;
    
     % put in periodic boundary
    MO(1:r,:)=MO(Order+1:Order+r,:);
    MO(Order+r+1:end,:)=MO(r+1:2*r,:);
    MO(:,1:r)=MO(:,Order+1:Order+r);
    MO(:,Order+r+1:end)=MO(:,r+1:2*r);
    
    N=sum(sum(MO(p(1)-r:p(1)+r,p(2)-r:p(2)+r)));
    NO=real(N);
    NG=imag(N);
    if og==1
        % Green itself
        NG=NG-1;
    else
        NO=NO-1;
    end     
end

function Plot(P,No,Order)
% Makes a plot of distribution of people
    plot(P(1:No,1),P(1:No,2),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
    hold on
    plot(P(No+1:end,1),P(No+1:end,2),'o','MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
    eps=1/(2*Order);
    X=(0:1/Order:1)-eps;
    for k=X
        plot([k,k],[-eps,1-eps],'k-')
    end
    for k=X
        plot([-eps,1-eps],[k,k],'k-')
    end
    axis equal
    set(gca,'visible','off')
    set(gcf,'color','w')
    hold off
end

function [Plist,Elist]=GenerateInitial7(InitialP1,InitialP2,InitialE,rr,re,No,Ng,Order,r,N)
    % Find Initial conditions 
    % Turn red and empty into blue ----Localized, random
    rep=0;
    maxrep=size(InitialE,3);
    Elist=zeros(size(InitialE,1)-re,2,maxrep);
    P1list=zeros(size(InitialP1,1)-rr,2,maxrep);
    P2list=zeros(size(InitialP2,1)+rr+re,2,maxrep);
    for t=1:size(InitialE,3)
        P1=InitialP1(:,:,t);
        P2=InitialP2(:,:,t);
        E=InitialE(:,:,t);
        % Find North/South people 
        P1N=P1(P1(:,2)>0.5,:);
        P1S=P1(P1(:,2)<=0.5,:);
        P2N=P2(P2(:,2)>0.5,:);
        P2S=P2(P2(:,2)<=0.5,:);
        EN=E(E(:,2)>0.5,:);
        ES=E(E(:,2)<=0.5,:);
        % Change North
        if size(EN,1)>re && size(P1N,1)>rr
            ib=rr+re;
            rep=rep+1;
            % randomly change red to empty
            EN=[EN;P1N(1:rr,:)];
            P1N=P1N(rr+1:end,:);
            EN=EN(randperm(size(EN,1)),:);
            % randomly introduce blue agents
            P2N=[P2N;EN(1:ib,:)];
            EN=EN(ib+1:end,:);
            % Put together
            P1F=[P1N;P1S];
            P2F=[P2N;P2S];
            EF=[EN;ES];
            P1list(:,:,rep)=P1F;
            P2list(:,:,rep)=P2F;
            Elist(:,:,rep)=EF;
        end
    end
    P1list=P1list(:,:,1:rep);
    P2list=P2list(:,:,1:rep);
    Elist=Elist(:,:,1:rep);
    tag=[zeros(No-rr,1);ones(Ng+rr+re,1)].*ones(1,1,rep);
    Plist=cat(2,cat(1,P1list,P2list),tag);
    
    % Equiliberate North
    tag=[zeros(No-rr,1);ones(Ng+rr+re,1)];
    parfor i=1:rep
        if mod(i,100)==0
            disp(i)
        end
        P=Plist(:,:,i);
        Empty=Elist(:,:,i);
        for j=1:2000
            List=1:size(P,1);
            List=List(P(:,2)>=0.5);
            pi=List(ceil(rand*length(List))); % random person north
            ListE=1:size(Empty,1);
            ListE=ListE(Empty(:,2)>=0.5);
            qi=N+ListE(ceil(rand*length(ListE))); % random empty spot
            % Find total utility/potential before switch
            U=FindUtility(pi,P,Order,r);
            Pi=Potential(pi,P);
            % Find P after switch
            LOC=[P(:,1:2);Empty];
            p=LOC(pi,:);
            q=LOC(qi,:);
            LOC(pi,:)=q;
            LOC(qi,:)=p;
            sP=[LOC(1:N,:),tag];
            sEmpty=LOC(N+1:end,:);
            V=FindUtility(pi,sP,Order,r);
            Pf=Potential(pi,sP);
            % Move with probability
            if rand<1/(1+exp((U+Pi)-(V+Pf)))
                P=sP;
                Empty=sEmpty;
            end
        end
        Plist(:,:,i)=P;
        Elist(:,:,i)=Empty;
    end
end
