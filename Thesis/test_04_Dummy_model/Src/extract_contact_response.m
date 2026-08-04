function y_contact = extract_contact_response(Struct, Interfaces, labels, q_constrained)
%EXTRACT_CONTACT_RESPONSE Storie temporali dei nodi di contatto, per interfaccia.
%
%   y_contact = EXTRACT_CONTACT_RESPONSE(Struct, Interfaces, labels, q)
%
%   Struct        oggetto AbaqusStructure (serve per unconstrain_vector)
%   Interfaces    struct dei metadati di interfaccia prodotta dal main
%   labels        cell array delle etichette da estrarre
%   q_constrained spostamenti sui GdL liberi, [n_dofs_c x n_time]
%
%   Restituisce una struct con un campo per interfaccia:
%       y_contact.<label>.X   [n_nodi x n_time]  spostamenti in X
%       y_contact.<label>.Y   [n_nodi x n_time]  spostamenti in Y
%
%   Entrambe le componenti vengono salvate a prescindere dalla direzione di
%   contatto dell'interfaccia: quella normale serve per il gap, quella
%   tangenziale per vedere lo scorrimento lungo la parete.

y_full    = Struct.AssemblyObj.unconstrain_vector(q_constrained);
y_contact = struct();

for i = 1:numel(labels)
    lbl = labels{i};
    if ~isfield(Interfaces, lbl) || isempty(Interfaces.(lbl).nodes)
        continue;
    end
    y_contact.(lbl).X = y_full(Interfaces.(lbl).global_X, :);
    y_contact.(lbl).Y = y_full(Interfaces.(lbl).global_Y, :);
end
end
