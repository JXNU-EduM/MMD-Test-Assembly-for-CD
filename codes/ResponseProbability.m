function result=ResponseProbability(model,ModelParameter,ItemMatrix,StudentMatrix)
    [ItemNumber,AttriNumber]=size(ItemMatrix);
    StudentNumber=size(StudentMatrix,1);
    outputMatrix=zeros(StudentNumber,ItemNumber);
    if strcmp('DINA',model)
        s=ModelParameter(1,:);
        g=ModelParameter(2,:);
       for i=1:StudentNumber
           for j=1:ItemNumber
               kx=1;
              for k=1:AttriNumber
                  kx=kx*StudentMatrix(i,k)^ItemMatrix(j,k);
              end 
              outputMatrix(i,j)=((1-s(j))^kx)*(g(j)^(1-kx));
           end
       end
       result=outputMatrix;
    elseif strcmp('DINO',model)
        s=ModelParameter(1,:);
        g=ModelParameter(2,:);
        for i=1:StudentNumber
            for j=1:ItemNumber
                omega=1;
                for k=1:AttriNumber
                    omega=omega*((1-StudentMatrix(i,k))^ItemMatrix(j,k));                
                end
                omega=1-omega;
                temp=((1-s(j))^omega)*(g(j)^(1-omega));
                outputMatrix(i,j)=temp;
            end
        end
        result=outputMatrix;
    elseif strcmp('R-RUM',model)
        pi=ModelParameter(1,:);
        ModelParameter(1,:)=[];
        r=ModelParameter';
        P=zeros(StudentNumber,ItemNumber);
        for i=1:StudentNumber
           for j=1:ItemNumber
               temp=1;
              for k=1:AttriNumber
                  temp=temp*r(j,k)^(ItemMatrix(j,k)*(1-StudentMatrix(i,k)));
              end
              P(i,j)=pi(1,j)*temp;
           end
        end
        result=P;
    else 
        result=-1;
    end
end
