function [result1,result2,result4]=SimDis(PrMatri,AttriMatrix)
    [StuNum,ItemNum]=size(PrMatri);
    [S,K]=size(AttriMatrix);
    Kpc=StuNum*(StuNum-1);
    SimpMatr=zeros(S);
    D=cell(1,ItemNum);
    CDI=zeros(1,ItemNum);
    h=zeros(StuNum);
    Vkl1=zeros(Kpc,ItemNum);
    Flag=ones(StuNum,1)*ones(1,StuNum)-eye(StuNum,StuNum);
    SKpc=K*2^(K-1);
    SVkl1=zeros(SKpc*2,ItemNum);
    for i=1:S
       for j=1:S
           temp1=AttriMatrix(i,:)==AttriMatrix(j,:);
           temp2=AttriMatrix(i,:)-AttriMatrix(j,:);
           if sum(temp1)==(K-1)&&sum(temp2)==1
              SimpMatr(i,j)=1; 
           end
       end
    end
    SimpMatr=SimpMatr+SimpMatr';
    for j=1:ItemNum
        D{j}=bsxfun(@times,log(PrMatri(:,j)*(1.0./PrMatri(:,j))'),PrMatri(:,j))+bsxfun(@times,log((1-PrMatri(:,j))*(1.0./(1-PrMatri(:,j))')),(1-PrMatri(:,j)));
        Vkl1(:,j)=D{j}(logical(Flag));
        SVkl1(:,j)=D{j}(logical(SimpMatr'));  
    end
    for u=1:StuNum
       for v=1:StuNum
           if u~=v
           for k=1:K
              h(u,v)=h(u,v)+(AttriMatrix(u,k)-AttriMatrix(v,k))^2;
           end
           h(u,v)=1.0/h(u,v);
           end
       end
    end
    for j=1:ItemNum
        Djuv=D{j};
        Dj=0;
        Dj=sum(sum(Djuv.*h));
        CDI(j)=Dj/sum(sum(h));
    end
    result1=CDI;
    result2=Vkl1;
    result4=SVkl1;
end






