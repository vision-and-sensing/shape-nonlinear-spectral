function [dTVmag_signed, dTVnormed_unsigned] = DeformationFieldSignedMagnitudesUnsignedDirections(F, V_orig, V_hat, do_draw, view_vec)
    % helper function
    rowwise_scalar_mult = @(u, v) sum(u .* v, 2);

    % displacement ("deformation") field
    dTV = V_orig - V_hat;

    % normals
    TR = triangulation(F, V_orig);
    N = vertexNormal(TR);

    % magnitude  of displacement field
    dTVmag = sqrt(dTV(:,1).^2 + dTV(:,2).^2 + dTV(:,3).^2);


    % normed displacement field
    dTVnormed = dTV ./ dTVmag;

    % sign of displacement field (indicates whether vectors point "in" or "out")
    dTVsign = sign(rowwise_scalar_mult(dTV, N));

    % unsigned normed displacement field
    dTVnormed_unsigned = dTVnormed .* dTVsign;

    % signed magnitude of displacement field
    dTVmag_signed = dTVmag .* dTVsign;
    
    % show
    if do_draw
        figure;
        trisurf(F,V_hat(:,1),V_hat(:,2),V_hat(:,3), dTVmag_signed);
        view ([view_vec]);
        shading interp;
        camlight left;
        axis equal;
        colormap jet
        colorbar;
        title('signed magnitude')
        drawnow
    end
end