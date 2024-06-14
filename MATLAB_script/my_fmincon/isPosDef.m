function tf = isPosDef(M)
    % return true if M is positive definite.
    % isSymm = issymmetric(M); % check whether H_k is symmetric or not
    d = eig(M); % calculate eigenvalues of H_k
    tf = all(d >= 0);
end