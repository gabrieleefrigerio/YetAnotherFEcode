function check_ortho(P, Mc, nm, nome)
    d  = sqrt(full(sum(P .* (Mc*P), 1)));    % lunghezze in norma massa
    Pn = P ./ d;                             % normalizzo: resta solo l'orientazione
    G  = full(Pn' * (Mc * Pn));  G = (G+G')/2;   % G(i,j) = cos(theta_ij)

    offd = abs(G - diag(diag(G)));
    [cmax, k] = max(offd(:));  [i,j] = ind2sub(size(G), k);

    fprintf('\n=== %s (r = %d) ===\n', nome, size(G,1));
    fprintf('spread lunghezze  : %.2e\n', max(d)/min(d));
    fprintf('coseno max        : %.4f  -> angolo %.2f deg  (coppia %d-%d)\n', ...
            cmax, acosd(min(cmax,1)), i, j);
    fprintf('coseno medio      : %.4f\n', mean(offd(~eye(size(G)))));
    fprintf('cond(G) collettivo: %.3e\n', cond(G));

    if nargin > 2 && ~isempty(nm)
        fprintf('  max cos dentro modali    : %.4f\n', max(max(offd(1:nm,1:nm))));
        fprintf('  max cos dentro statici   : %.4f\n', max(max(offd(nm+1:end,nm+1:end))));
        fprintf('  max cos modali-statici   : %.4f\n', max(max(offd(1:nm,nm+1:end))));
    end
end