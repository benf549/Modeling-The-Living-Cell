%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/14/21, coded on MATLAB _R2020b_ 
% 
%Takes in a randomly generated number between zero and 1 and a vector of
%frequencies. Computes the cumulative sum and steps through each index to
%identify the residue corresponds to the random number.
%
%Inputs: random number to convert to an index and list of frequencies corresponding to amino acid
%frequencies.
%Output: the index corresponding to the randomly generated number. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resi = getRes(num, frequencies)
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