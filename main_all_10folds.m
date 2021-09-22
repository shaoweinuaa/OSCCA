clear all
clc
close all
load KIRC.mat
currentFolder = pwd;
addpath(genpath(currentFolder))
indicator=data(:,2);
survivalTime=data(:,1)/max(data(:,1));
feature=data(:,4:168);
Image=feature(:,1:150); %150 image features
Gene=feature(:,151:165); %15 gene  features
opts.r1=1e-5;
opts.r2=1e-4;
opts.r3=1e-4;
sigma=1e-3;
split_number=4;
g_id=[];
i_id=[];
for i=1:10   % 10 fold feature selection
   i
  trainData=data(find(Indices~=i),:);
  [gene_id,image_id]=OSCCA(trainData,split_number,opts);
  g_id=[g_id;gene_id];
  i_id=[i_id;image_id];
end
    
AG=tabulate(g_id);
AI=tabulate(i_id);
g_index=find(AG(:,2)>=5);
i_index=find(AI(:,2)>=10);
saveFileData=[survivalTime,indicator,Gene(:,g_index),Image(:,i_index)];  
dlmwrite('selectedFeature.txt',saveFileData); 
       
  
       







