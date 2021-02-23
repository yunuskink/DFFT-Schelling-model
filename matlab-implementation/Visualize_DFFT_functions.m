% Takes the steady state counts of agents in each block from DataSheet and
% uses them to calculate the spatial vexation and global social frustration
% function

h=get_h(DataSheet,s); %Calculate joint histogram
[~,~,~,f,V1,V2,~,~,~,~,~,~,~,~] = TC_DFFT_Y_2(h,8,1,s); %Least squared fitting to extract parameters from histogram

figure;surf(f); colorbar;title("frustration function")
xlabel("# Blue agents");ylabel("# Red agents")
figure;imagesc(reshape(V1,[5,5]));colorbar;title("vexation for blue agents")
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
figure;imagesc(reshape(V2,[5,5]));colorbar; title("vexation for red agents")
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);

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
