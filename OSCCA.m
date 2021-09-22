function [gene_id,image_id]=OSCCA(data,split_number,opts)
currentFolder = pwd;
addpath(genpath(currentFolder))
r1=opts.r1;
r2=opts.r2;
r3=opts.r3;
n_H=150;
n_G=15;
indicator=data(:,2);
survivalTime=data(:,1);
feature=data(:,4:168);
feature=zscore(feature);
Image=feature(:,1:150);
Gene=feature(:,151:165);


XL=feature(find(indicator==1),:); %  Features for un-cenosored patients
YL=survivalTime(find(indicator==1)); %Survival Time for un-censored patients 

Image_L=XL(:,1:size(Image,2));
Gene_L=XL(:,size(Image,2)+1:size(Image,2)+size(Gene,2));

%Inequality Constrains for Genomic Data
[m_same_Gene,m_different_Gene]=divide_mean(survivalTime,Gene,indicator,split_number);
m_Gene=[m_same_Gene;m_different_Gene];   

%Inequality Constrains for Image Data
[m_same_Image,m_different_Image]=divide_mean(survivalTime,Image,indicator,split_number);
m_Image=[m_same_Image;m_different_Image];

%%%%%%%%%%%Initialize %%%%%%%%%%%%%%%%%%
w_H=zeros(n_H,1);
w_G=zeros(n_G,1);
for iter=1:10


J_G=zeros(n_G,1);
J_H=zeros(n_H,1);

Q_G=zeros(n_G,1);
Q_H=zeros(n_H,1);

R_G=zeros((split_number-1)*2,1);
R_H=zeros((split_number-1)*2,1);

rou1=1e-4;
rou2=1e-4;
theta_G=zeros((split_number-1)*2,1);
theta_H=zeros((split_number-1)*2,1);

A_H=2*r3*Image_L'*Image_L/size(XL,1);
B_H=-Image'*Gene*w_G/(size(Image,1)*size(Image,1)*10^2)+2*(r3/size(XL,1))*Image_L'*(Gene_L*w_G-YL);
C_H=m_Image;

%optimize W_H By ADMM
for i=1:500
J_H=lasso(eye(n_H),w_H+Q_H/rou1,'lambda',2*r1/rou1);
w_H=inv((A_H+rou1*eye(n_H)+rou1*C_H'*C_H))*(-B_H-Q_H+rou1*J_H-C_H'*R_H+rou1*C_H'*theta_H);
theta_H=max(zeros((split_number-1)*2,1),C_H*w_H+(1/rou1)*R_H);
Q_H=Q_H+rou1*(w_H-J_H);
R_H=R_H+rou2*(C_H*w_H-theta_H);
end

A_G=2*r3*Gene_L'*Gene_L/size(XL,1);
B_G=-Gene'*Image*w_H/((size(Gene,1)*size(Gene,1))*10^2)+2*(r3/size(XL,1))*Gene_L'*(Image_L*w_H-YL);
C_G=m_Gene; 

%optimize W_G By ADMM
for i=1:500
J_G=lasso(eye(n_G),w_G+Q_G/rou1,'lambda',2*r2/rou1);
w_G=inv((A_G+rou1*eye(n_G)+rou2*C_G'*C_G))*(-B_G-Q_G+rou1*J_G-C_G'*R_G+rou2*C_G'*theta_G);
theta_G=max(zeros((split_number-1)*2,1),C_G*w_G+(1/rou2)*R_G);
Q_G=Q_G+rou1*(w_G-J_G);
R_G=R_G+rou2*(C_G*w_G-theta_G);
end

WW_G(:,iter)=w_G;
WW_H(:,iter)=w_H;

end
gene_id=find(abs(w_G)>1e-4);
image_id=find(abs(w_H)>1e-4);

