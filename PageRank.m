function [pi, time, iter] = PageRank(h, n, alpha, epsilon)
    %INPUTS:
    %h = list of all hyperlinks in paired form
    %n = total number of sites
    %alpha = scaling parameter value
    %epsilon = convergence tolerance
    %OUTPUTS:
    %pi = PageRank vector
    %time = time taken to calculate PageRank vector
    %iter = number of iterations required for convergence
    %CODE:
    %h = sortrows(h);
    %split hyperlinks into lists of outlinks and inlinks
    %i = h(:,1);
    %j = h(:,2);
    %create sparse matrix of all hyperlinks. Storing in this form reduces memory allocation and complexity order of calculations. All values are equal to one, so the matrix is not yet stochastic.
    %H = sparse(i, j, 1, n, n);
    H=sparse(h);
    stochasticize = zeros(n, 1);
    a = zeros(n, 1);
    %for each hyperlink, increment outlinking page's link number by one
    for k = 1:length(i)
        stochasticize(i(k)) = stochasticize(i(k)) + 1;
    end
    counter = 0;
    for k = 1:n
        if stochasticize(k) ~= 0
            for m = 1:stochasticize(k)
                %divide link value by number of outlinks from corresponding page to make matrix stochastic
                counter = counter + 1;
                H(i(counter), j(counter)) = H(i(counter), j(counter)) / stochasticize(k);
            end
        else
            %add to dangling node vector if no outlinks from page
            a(k) = 1;
        end
    end
    %save on calculation time
    a = sparse(a);
    e = ones(n, 1);
    %set up initial PageRank vector
    pi = (1 / n) * transpose(e);
    iter = 0;
    %set up initial residual
    residual = 1;
    %start measuring time taken
    tic;

    while (residual >= epsilon)
        prevpi = pi;
        iter = iter + 1;
        %iterate PageRank equation on sparse vector; order o(nnz)
        pi = alpha * pi * H + (alpha * (pi * a) + 1 - alpha) * ((1 / n) * transpose(e));
        residual = norm(pi - prevpi, 1);
    end
    time = toc;
end