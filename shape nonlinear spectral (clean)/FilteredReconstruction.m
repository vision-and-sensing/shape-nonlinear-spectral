function [V_rec] = FilteredReconstruction(res, Phi, bands, amps)
    %% sizes
    assert(numel(bands) == numel(amps), 'assertion error, each band must recieve an amplification scalar')
    nv = size(res, 1);  % number of vertices
    nfn = size(res, 2);  % number of (scalar) functions consisting the vectorial function
    num_bands = numel(bands);
    %% sum phis in each band
    % summation helper function
    band_helper = @(band) squeeze(sum(reshape(transpose(cell2mat(Phi(band))), nv, nfn, 1+band(end)-band(1)), 3));    
    % sum each band
    bands_summed = zeros(nv, nfn, num_bands);
    for b = 1:num_bands
        bands_summed(:, :, b) = band_helper(bands{b});
    end

    %% sum amplified bands & reconstruct
    BPF = zeros(nv, nfn);
    for b = 1:num_bands
        BPF = BPF + amps(b) * bands_summed(:, :, b);
    end
    V_rec = res + BPF;
end