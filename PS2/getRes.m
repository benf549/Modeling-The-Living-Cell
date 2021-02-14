
function resi = getRes(num, frequencies)
%Takes in a randomly generated number between zero and 1 and a vector of
%frequencies. Computes the cumulative sum and steps through each index to
%identify the residue corresponds to the random number. 
    resi = 1;
    for i = cumsum(frequencies)
        if num < i
            break
        elseif num > 0.9989 %Because probabilities dont add to 1, add remaining ~0.001 to V
            resi = 20;
            break
        else
            resi = resi + 1;
        end
    end
end