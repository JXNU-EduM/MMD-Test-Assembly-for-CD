clear
clc
time=200;
Models={'DINA','DINO','R-RUM'};
RHOs=[0,0.5];
A=4;
J=20;
for A=4:4:8%attributes
    if A==4
        JM=[3 2 1;9 7 4];
    else
        JM=[4 3 2;9 7 4];
    end
    for m=1:1:3
        Model=Models{m};
        for r=1:2
            r
            rho=RHOs(r);  
            FunMain(time,Model,J,rho,A,JM,m)
        end
    end
end