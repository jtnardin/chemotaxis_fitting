function [A_pos,A_pos_0,A_pos_1,A_neg,A_neg_0,A_neg_1] = aMatrixConstruction(xn)

       %A matrices now must be space-dependent
    A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [vw.*(-1+sw/2); (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],xn,xn);
    
    A_pos_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*(1-1*se/2); ve.*se/2],xn,xn);

    A_pos_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
        (-vw.*sw/2)],total,total);



    A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [(-vw.*sw/2); (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],xn,xn);

    A_neg_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
        (vw.*sw/2-vw)],xn,xn);

    A_neg_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*se/2; (ve-ve.*se/2)],xn,xn);


end