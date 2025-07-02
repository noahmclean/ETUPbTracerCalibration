function F = vermeeshcopy3(S, n ,mu, smallS, datameas)
F = 0;
for i=1:n
    smallSi = smallS(:,:,i);
    A = inv([smallSi(1,1) + S(1,1), smallSi(1,2) + S(1,2), smallSi(1,3) + S(1,3)
             smallSi(2,1) + S(2,1), smallSi(2,2) + S(2,2), smallSi(2,3) + S(2,3)
             smallSi(3,1) + S(3,1), smallSi(3,2) + S(3,2), smallSi(3,3) + S(3,3)]);
    B = [datameas(1,i) - mu(1); datameas(2,i) - mu(2); datameas(3,i) - mu(3)];
    F = F + A*B*B'*A - A;
end