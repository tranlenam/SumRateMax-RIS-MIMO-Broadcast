function vect_out = pow_proj(Pt,vect_in)
    vect_in = vect_in.';
    [sort_val,sort_idx] = sort(vect_in,'descend');
    for n = length(vect_in):-1:1
        water_level = (sum(sort_val(1:n))-Pt)/n;
        di = sort_val(1:n)-water_level;
        if water_level<0
            disp('Negative water level');
        end
        if di(:)>=0
            break
        end   
    end
    vect_out = zeros(1,length(vect_in));
    vect_out(sort_idx(1:n)) = di;
end 