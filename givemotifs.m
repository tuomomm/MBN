%givemotifs(M): Calculate the number of each type of three-node motif.
%Tuomo M?ki-Marttunen, 2013-2016

function c=givemotifs(M)

N = size(M,1);
c = zeros(1,16); %number of found motifs 1-13
ctot = 0; %total number of found motifs



for v=1:N
    for v2=v+1:N
        for v3=v2+1:N
            n = M(v,v2) + M(v2,v) + M(v,v3) + M(v3,v) + M(v2,v3) + M(v3,v2); %number of connections in total in the triple
            if n==0                
                c(1) = c(1) + 1;
            elseif n==1
                c(2) = c(2) + 1;
            elseif n==2
                if M(v,v2) && M(v2,v) || M(v,v3) && M(v3,v) || M(v2,v3) && M(v3,v2)
                    c(4) = c(4) + 1;
                elseif M(v,v2) && M(v,v3) || M(v2,v) && M(v2,v3) || M(v3,v) && M(v3,v2)
                    c(6) = c(6) + 1;
                elseif M(v,v2) && M(v2,v3) || M(v,v3) && M(v3,v2) || M(v2,v3) && M(v3,v) || M(v2,v) && M(v,v3) || M(v3,v2) && M(v2,v) || M(v3,v) && M(v,v2)
                    c(5) = c(5) + 1;
                else %if M(v2,v) && M(v3,v) || M(v,v2) && M(v3,v2) || M(v,v3) && M(v2,v3)
                    c(3) = c(3) + 1;
                end
            elseif n==3
                if M(v,v2) && M(v2,v) && (M(v,v3) || M(v2,v3)) || ...
                   M(v,v3) && M(v3,v) && (M(v,v2) || M(v3,v2)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v2,v) || M(v3,v))
                    c(9) = c(9) + 1;
                elseif M(v,v2) && M(v2,v) && (M(v3,v) || M(v3,v2)) || ...
                   M(v,v3) && M(v3,v) && (M(v2,v) || M(v2,v3)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v,v2) || M(v,v3))
                    c(7) = c(7) + 1;
                elseif M(v,v2) && M(v2,v3) && M(v3,v) || M(v,v3) && M(v3,v2) && M(v2,v)
                    c(10) = c(10) + 1;
                else
                    c(8) = c(8) + 1;
                end
            elseif n==4
                if ~M(v,v2) && ~M(v2,v) || ~M(v,v3) && ~M(v3,v) || ~M(v2,v3) && ~M(v3,v2)
                    c(12) = c(12) + 1;
                elseif M(v,v2) && M(v2,v) && (M(v,v3) && M(v2,v3)) || ...
                   M(v,v3) && M(v3,v) && (M(v,v2) && M(v3,v2)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v2,v) && M(v3,v))
                    c(14) = c(14) + 1;
                elseif M(v,v2) && M(v2,v) && (M(v3,v) && M(v3,v2)) || ...
                   M(v,v3) && M(v3,v) && (M(v2,v) && M(v2,v3)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v,v2) && M(v,v3))
                    c(11) = c(11) + 1;
                else
                    c(13) = c(13) + 1;
                end
            elseif n==5
                c(15) = c(15) + 1;
            else
                c(16) = c(16) + 1;
            end
        end
    end
end
