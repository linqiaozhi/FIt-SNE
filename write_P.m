
function write_P(P)
% WRITE_P Writes a sparse matrix P to files that FIt-SNE can read in lieu of computing P based on nearest neighbors. The three vectors it produces are as follows: 
% val_P:  The K non-zero values of a sparse N by N matrix P
% col_P: a K-length vector of unsigned ints, giving the column indices of each element
% row_P: a vector of N+1 length of unsigned ints, giving the index in the preceding two matrices corresponding to where  each row "starts".

    delete temp/*.dat
    [N,~] = size(P);
    row_P = zeros(N+1,1);
    row_P(2:end) = full(cumsum(sum(P~=0,2)));
    [I, J, V] = find(P');
    val_P = V;
    col_P = I-1;
    h = fopen('temp/val_P.dat', 'wb');
    fwrite(h, val_P, 'double');
    fclose(h);
    h = fopen('temp/col_P.dat', 'wb');
    fwrite(h, col_P, 'integer*4');
    fclose(h);
    h = fopen('temp/row_P.dat', 'wb');
    fwrite(h, row_P, 'integer*4');
    fclose(h);
end
