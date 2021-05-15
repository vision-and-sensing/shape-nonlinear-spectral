function [V_rec] = FilteredReconstructionMatPhi(res, phi, fltr)
    nfn = size(res, 2);  % number of (scalar) functions consisting the vectorial function
    BPF = zeros(size(res));
    for f = 1:nfn
        BPF(:, f) = sum(squeeze(phi(:, f, :)) .* fltr, 2);
    end
    V_rec = res + BPF;
end