function padj = palm_fastfdr(pval)

% Sort p-values
V = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR] = sort(oidx);

% Loop over each sorted p-value
padj = zeros(size(pval));
prev = 1;
for i = V:-1:1,
    % The p-adjusted for the current p-value is the smallest slope among
    % all the slopes of each of the p-values larger than the current one
    % Yekutieli & Benjamini (1999), equation #3. Note that here c(V) = 1.
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);