function FunMain(time,Mo,J,rho,AttrNum,JMatrix,Moindex)
    stdy=0;
    StuNum=10000;
    ItemBNum=300;
    RaCDI=zeros(6,time);
    RaKl=zeros(6,time);
    RaSKl=zeros(6,time);
    RaDiff=zeros(6,time);
    RaSDiff=zeros(6,time);
    RaRand=zeros(6,time);
    TimeCount=zeros(4,time);
    for i=1:time
        [StuMatr1,ItemMatrix1,AttrMatr1,StuMatrIndex1]=SimulData(AttrNum,ItemBNum,StuNum,rho,stdy);
        [StuMatr,ItemMatrix,AttrMatr,StuMatrIndex]=SimulData(AttrNum,ItemBNum,StuNum,rho,stdy);
        ParMatr=ModelParameter(Mo,AttrNum,ItemBNum);
        if strcmp('R-RUM',Mo)
            ParMatr(2:AttrNum+1,:)=ParMatr(2:AttrNum+1,:).*ItemMatrix';
        end
        PMatrix=ResponseProbability(Mo,ParMatr,ItemMatrix,AttrMatr);
        U=rand(StuNum,ItemBNum)<PMatrix(StuMatrIndex,:);    
        prior=mean(bsxfun(@eq,StuMatrIndex1,1:2^AttrNum))';
        [Dj,Dkl,SDkl]=SimDis(PMatrix,AttrMatr);
        [Rm]=RandCon(Dj);
        [NC,IC,AC,LNC,LIC,LAC,TC]=CDI(Dj,ItemMatrix,J,JMatrix);
        [RaCDI(1,i),RaCDI(4,i)]=map(U(:,LNC),PMatrix(:,LNC),prior,AttrMatr,StuMatr,StuMatrIndex);
        [RaCDI(2,i),RaCDI(5,i)]=map(U(:,LIC),PMatrix(:,LIC),prior,AttrMatr,StuMatr,StuMatrIndex);
        [RaCDI(3,i),RaCDI(6,i)]=map(U(:,LAC),PMatrix(:,LAC),prior,AttrMatr,StuMatr,StuMatrIndex);
        tic
        [SNK,SIK,SAK,SLNK,SLIK,SLAK,TSK]=MiLp(SDkl,ItemMatrix,J,JMatrix);  
        [RaSKl(1,i),RaSKl(4,i)]=map(U(:,SLNK),PMatrix(:,SLNK),prior,AttrMatr,StuMatr,StuMatrIndex);
        [RaSKl(2,i),RaSKl(5,i)]=map(U(:,SLIK),PMatrix(:,SLIK),prior,AttrMatr,StuMatr,StuMatrIndex);
        [RaSKl(3,i),RaSKl(6,i)]=map(U(:,SLAK),PMatrix(:,SLAK),prior,AttrMatr,StuMatr,StuMatrIndex);
        TimeCount(1,i)=toc;
        tic
        if(AttrNum==4)
            [SNK,SIK,SAK,LNK,LIK,LAK,TK]=MiLp(Dkl,ItemMatrix,J,JMatrix);
            [RaKl(1,i),RaKl(4,i)]=map(U(:,LNK),PMatrix(:,LNK),prior,AttrMatr,StuMatr,StuMatrIndex);
            [RaKl(2,i),RaKl(5,i)]=map(U(:,LIK),PMatrix(:,LIK),prior,AttrMatr,StuMatr,StuMatrIndex);
            [RaKl(3,i),RaKl(6,i)]=map(U(:,LAK),PMatrix(:,LAK),prior,AttrMatr,StuMatr,StuMatrIndex);
        end
        TimeCount(2,i)=toc;
        [NR,IR,AR,LNR,LIR,LAR,TR]=MiLp(Rm,ItemMatrix,J,JMatrix);
         [RaRand(1,i),RaRand(4,i)]=map(U(:,LNR),PMatrix(:,LNR),prior,AttrMatr,StuMatr,StuMatrIndex);
         [RaRand(2,i),RaRand(5,i)]=map(U(:,LIR),PMatrix(:,LIR),prior,AttrMatr,StuMatr,StuMatrIndex);
         [RaRand(3,i),RaRand(6,i)]=map(U(:,LAR),PMatrix(:,LAR),prior,AttrMatr,StuMatr,StuMatrIndex);
         if (AttrNum==4)
            Cri=[AttrNum*ones(6,1) Moindex*ones(6,1) rho*ones(6,1)  i*ones(6,1) [1*ones(3,1);2*ones(3,1)] [1:3 1:3]' RaRand(:,i) RaCDI(:,i) RaKl(:,i)  RaSKl(:,i) TR,  TC,TK,TSK];
         else
            Cri=[AttrNum*ones(6,1) Moindex*ones(6,1) rho*ones(6,1)  i*ones(6,1) [1*ones(3,1);2*ones(3,1)] [1:3 1:3]' RaRand(:,i) RaCDI(:,i) nan*ones(6,1)  RaSKl(:,i) TR,  TC,nan*ones(6,1),TSK];
         end
         dlmwrite('study2_200_8_T_20201205.txt',Cri,'-append');
    end
end