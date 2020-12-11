function [result1,result2,result3,result4]=SimulData(An,In,Sn,rho,StDy)
    mu=zeros(An,1);
    sigma=eye(An);
    for i=1:An
        for j=1:An
            if i~=j
                sigma(i,j)=rho;
            end
        end
    end
    r1=mvnrnd(mu,sigma,Sn);
    for i=1:Sn
        for j=1:An
            if StDy==0
                if r1(i,j)>=0
                    r1(i,j)=1;
                else
                    r1(i,j)=0;
                end
            else
                Ci=normcdf(i/(i+1),0,1);
                if r1(i,j)>=Ci
                    r1(i,j)=1;
                else
                    r1(i,j)=0;
                end
            end
        end
    end
    dim=An;
    inputMatrix=dec2bin(0:2^dim-1);
    outputMatrix=zeros(2^dim,dim);
    for i=1:2^dim
        temp=inputMatrix(i,:);
        outputMatrix(i,:)=str2num(temp(:));
    end
    outputMatrix(1,:)=[];
    tempRand=randi(2^An-1,1,In);
    ItemBankMatrix=zeros(In,An);
    for i=1:In
        ItemBankMatrix(i,:)=outputMatrix(tempRand(i),:);
    end
    r2=ItemBankMatrix;
    dim=An;
    r3intput=dec2bin(0:2^dim-1);
    r3=zeros(2^dim,dim);
    for i=1:2^dim
        temp=r3intput(i,:);
        r3(i,:)=str2num(temp(:));
    end
    result1=r1;
    result2=r2;
    result3=r3;
    result4=zeros(Sn,1);
    for i=1:Sn
        for k=1:2^dim
            if isequal(result1(i,:),result3(k,:))
                result4(i,1)=k;
            end
        end
    end
end