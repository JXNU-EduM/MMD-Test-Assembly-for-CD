function [result1,result2]=SimTest(model,ParaMatrix,Item,SMatrix,TestAss,Qs)
    [ItemNumber,AttriNumber]=size(Item); 
    [StuNumber,~]=size(SMatrix);
    estS=ones(StuNumber,2^AttriNumber);
    estStu=zeros(StuNumber,AttriNumber);
    if strcmp('DINA',model)||strcmp('DINO',model)
        TestPara=zeros(2,ItemNumber);
        for i=1:ItemNumber
            TestPara(:,i)=ParaMatrix(:,TestAss(i));
        end
        zerosPara=zeros(2,ItemNumber);
        RespMatrix=ResponseProbability(model,zerosPara,Item,SMatrix);
        SimPesProb=ResponseProbability(model,TestPara,Item,SMatrix);
        Pj=ResponseProbability(model,TestPara,Item,Qs);
    elseif strcmp('R-RUM',model)
        TestPara=zeros(AttriNumber+1,ItemNumber);
        for i=1:ItemNumber
           TestPara(:,i)=ParaMatrix(:,TestAss(i));
        end
        SimPesProb=ResponseProbability(model,TestPara,Item,SMatrix);
        Pj=ResponseProbability(model,TestPara,Item,Qs);
    end
    for i=1:StuNumber
       for j=1:ItemNumber
           tr=rand(1);
          if SimPesProb(i,j)>=tr;
              RespMatrix(i,j)=1;
          else
              RespMatrix(i,j)=0;         
          end
       end
    end
    for i=1:StuNumber
       for j=1:ItemNumber
           for k=1:2^AttriNumber
              temp1=Pj(k,j)^RespMatrix(i,j);
              temp2=(1-Pj(k,j))^(1-RespMatrix(i,j));
              estS(i,k)=estS(i,k)*temp1*temp2; 
           end
       end  
        [maxest,index]=max(estS(i,:));
        estStu(i,:)=Qs(index,:);
    end
    Ra=0;
    for i=1:StuNumber
        L=estStu(i,:)==SMatrix(i,:);
        if sum(L)==AttriNumber
            Ra=Ra+1;
        end
    end
    Rap=Ra/StuNumber;
    result1=mean(mean( estStu==SMatrix));
    result2=Rap;
end