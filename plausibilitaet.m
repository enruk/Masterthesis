function [wrong_seq,wrong_allo] = plausibilitaet (Population,p,n)

wrong_seq = zeros(1,n);
wrong_allo = zeros(1,n);

index = 1;
for z1=1:p
    plausibel_seq = unique(Population{z1,1}(1,:));
    if length(plausibel_seq) ~= n
        wrong_seq(1,index) = z1;
        index = index + 1;
    end
end

wrong_seq(wrong_seq==0)=[];


index = 1;
for z1=1:p
    for z2=1:n
        if Population{z1,1}(2,z2) == Population{z1,1}(3,z2)
            wrong_allo(1,index) = z1;
            index = index + 1;
        end
    end
end

wrong_allo(wrong_allo==0)=[];


