function Bk = force_symmetry(Bk0)
for iaa=1:length(Bk0)
    for ibb=1:iaa
        Bk(iaa,ibb)=Bk0(iaa,ibb);
        Bk(ibb,iaa)=Bk0(iaa,ibb);
%         (yk*yk')/(yk'*sk)-(Bk*sk*sk'*Bk')/ (sk'*Bk*sk)
    end
end
end