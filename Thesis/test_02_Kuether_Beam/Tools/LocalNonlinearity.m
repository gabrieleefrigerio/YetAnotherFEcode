classdef LocalNonlinearity < handle
    % LOCALNONLINEARITY Models the unilateral spring
    
    properties
        k       % Spring stiffness 
        a       % Clearance 
        dof_idx % Index of the contact DOF
    end
    
    methods
        function obj = LocalNonlinearity(stiffness, clearance, dof_index)
            obj.k = stiffness;
            obj.a = clearance;
            obj.dof_idx = dof_index;
        end
        
        function f_nl = evaluate(obj, x, num_total_dofs)
            % Evaluates the non-linear force vector f_NL(x) in the full space
            f_nl = zeros(num_total_dofs, 1);
            displacement = x(obj.dof_idx);
            
            % Piecewise-linear logic (if it hits the spring, the force acts)
            if displacement > obj.a
                f_nl(obj.dof_idx) = obj.k * (displacement - obj.a);
            end
        end
        
        function f_nl_red = evaluate_reduced(obj, q, P)
            % f_NL projected into the reduced space: P' * f_NL(P*q)
            x_full = P * q;
            f_nl_full = obj.evaluate(x_full, size(P,1));
            f_nl_red = P' * f_nl_full;
        end
        
        function J_nl_red = jacobian_reduced(obj, q, P)
            % Computes the Jacobian projected into the reduced space: P' * J_NL * P
            x_full = P * q;
            displacement = x_full(obj.dof_idx);
            
            num_total_dofs = size(P, 1);
            J_nl_full = spalloc(num_total_dofs, num_total_dofs, 1); % Empty sparse matrix
            
            % If the spring is active, the tangent stiffness is k
            if displacement > obj.a
                J_nl_full(obj.dof_idx, obj.dof_idx) = obj.k;
            end
            
            % Projection
            J_nl_red = P' * J_nl_full * P;
        end
    end
end