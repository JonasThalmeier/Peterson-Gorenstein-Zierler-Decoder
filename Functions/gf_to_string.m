function [out_str] = gf_to_string(gf_in,alpha,m,mode)
%GF_TO_STRING Summary of this function goes here
%   Detailed explanation goes here
N = length(gf_in.x);
out_str = "";
if mode == 1
    for idx = 1:N
        if gf_in.x(idx) == 0
            continue;
        end
        exp = log(gf(gf_in.x(idx),m))./log(alpha);
        out_str = strcat(out_str, "α^", num2str(exp), "x^", num2str(N-idx));
        if idx ~= N
            out_str = strcat(out_str,'+');
        end
    end
else if mode == 2
        for idx = 1:N
            if gf_in.x(idx) == 0
                char = "0";
            elseif gf_in.x(idx) == alpha.x
                char = "α";
            elseif gf_in.x(idx) == 1
                char = "1";
            else
                exp = log(gf(gf_in.x(idx),m))./log(alpha);
                char = strcat("α^", num2str(exp));
            end
            out_str = strcat(out_str, char);
            if idx ~= N
                out_str = strcat(out_str,',');
            end
        end
else
    for idx = 1:N
        out_str = strcat(out_str, dec2bin(gf_in.x(idx),m));
        if idx ~= N
            out_str = strcat(out_str,',');
        end
    end
end
end
