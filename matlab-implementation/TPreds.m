% Calculates OMove for prediction timescale

% OUTPUT
% OMove: The number of moves in the underlying DFFT model for the same 
% number of trials as Schelling1
%% Extract frustration (calculate derivative) and vexation
h=get_h(DataSheet,s);
[~,~,~,f,V1,V2,~,~,~,~,~,~,~,~] = TC_DFFT_Y_2(h,8,1,s);
gf=NaN*zeros(s+1);
gf(1:size(f,1),1:size(f,2))=f;
f=gf;
[D1,D11,D2,D22]=Derivative2(f);
PMove=0;

%% Start simulation
M=iNo;
N=iNg;
E=Order^2-M-N;

i=0;
n1=zeros(Nblock,1);
n2=n1;
for m=1:CellSize % Row
    for n=1:CellSize % Column
        i=i+1;
        % Find number of orange and green people:
        n1(i)=DataSheet{m,n}(end,1);
        n2(i)=DataSheet{m,n}(end,2);
    end
end

% Initial dist of agent (list of blocks for each agent):
P1=[];
for i=1:Nblock
    app=i*ones(n1(i),1);
    P1=[P1;app];
end

P2=[];
for i=1:Nblock
    app=i*ones(n2(i),1);
    P2=[P2;app];
end

% initial list of empty cell (list of blocks for each cell)
e=s-n1-n2;
Ept=[];
for i=1:Nblock
    app=i*ones(e(i),1);
    Ept=[Ept;app];
end
    
for i=1:SampleSize
    if mod(i,1000)==0
        disp(i)
    end
    % generate list of agents for one cycle to move from
    List=ceil((M+N)*rand(SampleStep,1))';
    % generate list of empty cells to move to
    List2=ceil(E*rand(SampleStep,1))';
    
    for j=1:SampleStep
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
            PMove=PMove+1;
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
end
    
    

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