function [m_same,m_different]= divide_mean(survivalTime,feature,indicator,splitNumber)
index_death=find(indicator==1);
m_live=zeros(splitNumber,size(feature,2));
m_dead=zeros(splitNumber,size(feature,2));
[A,sortIndex]=sort(survivalTime);
pieceSize=ceil(length(survivalTime)/splitNumber);
tempIndex=0;
for i=1:splitNumber
 if i~=splitNumber
   intervalIndex=sortIndex(tempIndex+1:tempIndex+pieceSize);
   tempIndex=tempIndex+pieceSize;
 else
   intervalIndex=sortIndex(tempIndex+1:length(indicator));  
 end
   intervalIndexDeath=intersect(intervalIndex,index_death);
   intervalIndexLive=setdiff(intervalIndex,intervalIndexDeath);
   m_live(i,:)=mean(feature(intervalIndexLive,:));
   m_dead(i,:)=mean(feature(intervalIndexDeath,:)); 
end

for i=2:splitNumber
   m_same(i-1,:)=m_dead(i,:)-m_dead(i-1,:);
   m_different(i-1,:)=m_live(i,:)-m_dead(i-1,:);
end