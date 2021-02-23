% Evolve the altered ensemble given by Plist and Elist in Schelling1;
% Convert Plist and Elist into P1list and P2list to be used in prediction

% OUTPUTS:
% Data1: Time evolution of ensemble (red). Data(i,j,k) gives the number of 
% red agents at time t (frame) for j-th block in k-th city in the ensemble  
% Data2: Time evolution of ensemble (blue). Data(i,j,k) gives the number of 
% blue agents at time t (frame) for j-th block in k-th city in the ensemble 

% P1list(:,i) is the list of block indices of all red agents for the i-th 
% city in the altered ensemble (the same ensemble as Plist)
% P2list(:,i) is the list of block indices of all red agents for the i-th 
% city in the altered ensemble (the same ensemble as Plist)
%% Parameters:

T=150; % Time frames recorded
step=1000; % Steps between time samples

%% Specify Potential that matches potential.m (for faster speed)
[y,x]=meshgrid(0:1/Order:1-1/Order);
PO=-1.3*x; % Red spatial utility
PG=-1.3*y; % Blue spatial utility


%% Main
EL=Order/CellSize;
Data1=zeros(T,Nblock,SampleSize);
Data2=Data1;
No=iNo-rr;
Ng=iNg+rr+re; 
N=No+Ng;
ep=1/(2*Order); % correct for float point error when comparing sizes

% Calculate Cell Neighbor Index
[idx]=FindNBindex(Order);


% Initialize Various variables
BlkIdx=reshape(1:CellSize^2,CellSize,CellSize);
BlkIdx=BlkIdx(:,end:-1:1);
DataM=zeros(CellSize,Order);
DataMf=zeros(CellSize);


for rep=1:SampleSize
    if mod(rep,100)==0
        disp(rep)
    end
    
    % initialize positions: Each row contains to x,y position of one person
    P=Plist(:,:,rep);
    Empty=Elist(:,:,rep);
    LEmpty=length(Empty);
    tag=[zeros(No,1);ones(Ng,1)];
    
    % Convert into matrix representation (inverted matrix from actual coordinates)
    MO=zeros(Order);
    MG=zeros(Order);
    O=P(P(:,3)==0,1:2);
    O=round(O*Order)+1;
    Oindex=O(:,1)+(O(:,2)-1)*Order;
    MO(Oindex)=1;
    G=P(P(:,3)==1,1:2);
    G=round(G*Order)+1;
    Gindex=G(:,1)+(G(:,2)-1)*Order;
    MG(Gindex)=1;
    
    % Find list representation
    LP=[[Oindex;Gindex],tag];
    E=round(Empty*Order)+1;
    LE=E(:,1)+(E(:,2)-1)*Order;

    % Count # of agents in each block from their coordinates 
    for i=1:CellSize
        DataM(i,:)=sum(MO((i-1)*EL+1:i*EL,:),1);
    end
    
    for i=1:CellSize
        DataMf(:,i)=sum(DataM(:,(i-1)*EL+1:i*EL),2);
    end
    Data1(1,:,rep)=reshape(DataMf,1,[]);
    
    for i=1:CellSize
        DataM(i,:)=sum(MG((i-1)*EL+1:i*EL,:),1);
    end
    
    for i=1:CellSize
        DataMf(:,i)=sum(DataM(:,(i-1)*EL+1:i*EL),2);
    end
    Data2(1,:,rep)=reshape(DataMf,1,[]);

    % Run
    for j=2:T
        for k=1:step
            % Randomly pick a pair of indices 
            pi=ceil(rand*N); % random person
            qi=ceil(rand*LEmpty); % random empty spot
            og=LP(pi,2);
            pidx=LP(pi,1);
            qidx=LE(qi);
            % Find total utility/potential before switch
            U=FindUtility(pidx,og,MO,MG,idx);
            Pi=FindPotential(pidx,og,PO,PG);
            % Find MO and MG after switch
            MOf=MO;
            MGf=MG;
            if og
                % Green Move
                MGf(qidx)=1;
                MGf(pidx)=0;
            else
                MOf(qidx)=1;
                MOf(pidx)=0;
            end
            V=FindUtility(qidx,og,MOf,MGf,idx);
            Pf=FindPotential(qidx,og,PO,PG);
            % Move with probability
            if rand<1/(1+exp((U+Pi)-(V+Pf)))
                LP(pi,1)=qidx;
                LE(qi)=pidx;
                MO=MOf;
                MG=MGf;
            end
        end
        
        % Sort result into DataSheet
        % Count # of agents in each block from their coordinates 
        for i=1:CellSize
            DataM(i,:)=sum(MO((i-1)*EL+1:i*EL,:),1);
        end
        
        for i=1:CellSize
            DataMf(:,i)=sum(DataM(:,(i-1)*EL+1:i*EL),2);
        end
        Data1(j,:,rep)=reshape(DataMf,1,[]);
        
        for i=1:CellSize
            DataM(i,:)=sum(MG((i-1)*EL+1:i*EL,:),1);
        end
        
        for i=1:CellSize
            DataMf(:,i)=sum(DataM(:,(i-1)*EL+1:i*EL),2);
        end
        Data2(j,:,rep)=reshape(DataMf,1,[]);
    end
end

Data1=Data1(:,reshape(BlkIdx,1,[]),:);
Data2=Data2(:,reshape(BlkIdx,1,[]),:);

%% Create altered ensemble for use in predictions 
[P1list,P2list]=InitialCondSim(Data1,Data2); 

%% Plot Observed time evolution
Scale=1;
tf=T; % Change the plotted number of frames here (<=T)
Blockid=25;  % Change block number here
Pt=zeros(s+1,s+1,tf);
mean1=zeros(tf,1);
mean2=mean1;
for t=1:tf
    Data=[reshape(Data1(t,Blockid,:),[],1,1),...
        reshape(Data2(t,Blockid,:),[],1,1)];
    H=hist3(Data,'Ctrs',{0:s 0:s},'CdataMode','auto');
    P=H/sum(sum(H));
    Pt(:,:,t)=P;
    mean1(t)=nanmean(Data(:,1));
    mean2(t)=nanmean(Data(:,2));
end

% Find 1d pd
P1=reshape(sum(Pt,2),s+1,[],1);
P2=reshape(sum(Pt,1),s+1,[],1);
figure
imagesc(Scale*(0:tf-1),0:s,P1)
ax=gca;
ax.YDir='normal';
colorbar
colormap bone
hold on
plot(Scale*(0:tf-1),mean1,'r-','LineWidth',2)
figure
imagesc(Scale*(0:tf-1),0:s,P2)
ax=gca;
ax.YDir='normal';
colormap bone
hold on
plot(Scale*(0:tf-1),mean2,'r-','LineWidth',2)

%% Functions 
function [u]=FindUtility(p,og,MO,MG,idx)
% Find utility of agent with index p
    Neighbors=idx(idx(:,1)==p,2:end);
    NO=sum(MO(Neighbors));
    NG=sum(MG(Neighbors));
    u=Utility(NO,NG,og);
end

function [v]=FindPotential(p,og,PO,PG)
    if og
        % green
        v=PG(p);
    else
        v=PO(p);
    end
end

function [idx]=FindNBindex(CellSize) 
    % Finds indices for neighboring blocks/cells
    M=reshape(1:CellSize^2,CellSize,CellSize);
    ME=[M(end,:);M;M(1,:)];
    ME=[ME(:,end),ME,ME(:,1)];
    c1=reshape(ME(1:CellSize,1:CellSize)',CellSize^2,1);
    c2=reshape(ME(1:CellSize,2:CellSize+1)',CellSize^2,1);
    c3=reshape(ME(1:CellSize,3:CellSize+2)',CellSize^2,1);
    c4=reshape(ME(2:CellSize+1,1:CellSize)',CellSize^2,1);
    c5=reshape(ME(2:CellSize+1,2:CellSize+1)',CellSize^2,1);
    c6=reshape(ME(2:CellSize+1,3:CellSize+2)',CellSize^2,1);
    c7=reshape(ME(3:CellSize+2,1:CellSize)',CellSize^2,1);
    c8=reshape(ME(3:CellSize+2,2:CellSize+1)',CellSize^2,1);
    c9=reshape(ME(3:CellSize+2,3:CellSize+2)',CellSize^2,1);
    idx=[c5,c1,c3,c7,c9,c2,c4,c6,c8];
end

function [P1list,P2list]=InitialCondSim(Data1,Data2)
    % Obtain initial conditions for the Barker simulation


    P1list=zeros(sum(Data1(1,:,1)),size(Data1,3));
    for rep=1:size(Data1,3)
        list=[];
        for i=1:size(Data1,2)
            app=i*ones(Data1(1,i,rep),1);
            list=[list;app];
        end
        P1list(:,rep)=list;
    end

    P2list=zeros(sum(Data2(1,:,1)),size(Data2,3));
    for rep=1:size(Data2,3)
        list=[];
        for i=1:size(Data2,2)
            app=i*ones(Data2(1,i,rep),1);
            list=[list;app];
        end
        P2list(:,rep)=list;
    end
end

