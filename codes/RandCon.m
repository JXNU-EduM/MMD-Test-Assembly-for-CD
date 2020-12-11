function result=RandCon(Dis)
    [RowDis,ColumnDis]=size(Dis);
    MaxDis=max(Dis);
    MinDis=min(Dis);
    out=MinDis+(MaxDis-MinDis)*rand(RowDis,ColumnDis);
    result=out;
end