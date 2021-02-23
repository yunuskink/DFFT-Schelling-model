% Make Unscaled predictions of time evolution

% OUTPUT
% Data1 and Data2 in the same format as in Schelling2. 

%% Extract frustration (calculate derivative) and vexation
h=get_h(DataSheet,s); %Calculate joint histogram
[~,~,~,f,V1,V2,~,~,~,~,~,~,~,~] = TC_DFFT_Y_2(h,8,1,s); %Least squared fitting to extract parameters from histogram
gf=NaN*zeros(s+1);
gf(1:size(f,1),1:size(f,2))=f;
f=gf;
[D1,D11,D2,D22]=Derivative2(f);


%% Start simulation
M=iNo-rr;
N=iNg+rr+re;
E=Order^2-M-N;
Data1=zeros(T,Nblock,size(P1list,2));
Data2=Data1;

for rep=1:size(P1list,2)
    if mod(rep,100)==0
        disp(rep)
    end
    % Initial dist of agent (list of blocks for each agent):
    P1=P1list(:,rep);
    P2=P2list(:,rep);
    
    % initial agent count in each block
    [n1,~]=histcounts(P1,1:Nblock+1); 
    [n2,~]=histcounts(P2,1:Nblock+1);
    
    % initial list of empty cell (list of blocks for each cell)
    e=s-n1-n2;
    Ept=[];
    for i=1:Nblock
        app=i*ones(e(i),1);
        Ept=[Ept;app];
    end
    
    Data1(1,:,rep)=n1;
    Data2(1,:,rep)=n2;

    for i=2:T
        % generate list of agents for one cycle to move from
        List=ceil((M+N)*rand(step,1))';
        % generate list of empty cells to move to
        List2=ceil(E*rand(step,1))';
        
        for j=1:step
            k=List(j);
            q=Ept(List2(j));
            if k<=M
                %Red
                p=P1(k);
                potp=D1(n1(p)+1,n2(p)+1)+V1(p);
                potq=D11(n1(q)+1,n2(q)+1)+V1(q);
            else
                %Blue
                p=P2(k-M);
                potp=D2(n1(p)+1,n2(p)+1)+V2(p);
                potq=D22(n1(q)+1,n2(q)+1)+V2(q);
            end
            if (1/(exp(potq-potp)+1))>rand
                Ept(List2(j))=p;
                %move
                if k<=M
                    %red
                    P1(k)=q;
                    n1(q)=n1(q)+1;
                    n1(p)=n1(p)-1;
                else
                    %blue
                    P2(k-M)=q;
                    n2(q)=n2(q)+1;
                    n2(p)=n2(p)-1;
                end
            end
        end        
        Data1(i,:,rep)=n1;
        Data2(i,:,rep)=n2;
    end
end

%% Plot Predicted time evolution (With time scaling)

% Scale=PMove/OMove;
Scale = 1;
tf=T;
Blockid=25; % Change block number here
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
plot(Scale*(0:tf-1),Mean1(1:tf,Blockid),'b:','LineWidth',2)
figure
imagesc(Scale*(0:tf-1),0:s,P2)
ax=gca;
ax.YDir='normal';
colorbar
colormap bone
hold on
plot(Scale*(0:tf-1),mean2,'r-','LineWidth',2)
plot(Scale*(0:tf-1),Mean2(1:tf,Blockid),'b:','LineWidth',2)


%% Functions
function [D1,D11,D2,D22]=Derivative2(f)
    % Take partial derivatives (scaled) of f along first and second dimension
    F=[nan(1,size(f,2));f;nan(1,size(f,2))];
    ad=[f;nan(2,size(f,2))];
    subt=[nan(2,size(f,2));f];
    D11=ad-F;
    D1=F-subt;
    D1=D1(2:size(D1,1)-1,:);
    D11=D11(2:size(D11,1)-1,:);
    
    F=[nan(size(f,2),1),f,nan(size(f,2),1)];
    ad=[f,nan(size(f,1),2)];
    subt=[nan(size(f,1),2),f];
    D22=ad-F;
    D2=F-subt;
    D2=D2(:,2:size(D2,2)-1);
    D22=D22(:,2:size(D22,2)-1);
end

function [h] = get_h(DataSheet,s)
%Get the joint histograms
    CellSize=size(DataSheet,1);
    h=zeros(s+1,s+1,CellSize^2); 
    k=0;
    for i=1:CellSize %for each bin
        for j=1:CellSize
            k=k+1;
            Data=DataSheet{i,j};
            h(:,:,k)=hist3(Data,'Ctrs',{0:s 0:s},'CdataMode','auto');
        end
    end
end