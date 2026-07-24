function varargout=vout(C)
% small utility function for calculations on S.IE fields without endloos
% loops over all internal elements
% https://nl.mathworks.com/matlabcentral/answers/67874-generate-a-comma-separated-list-from-a-numeric-array
% trick is to make a vector of all values of say IE(:).r , do the
% calculations, and distribute again by this function. 

    if isnumeric(C), C=num2cell(C); end
    C=C(:).';
    varargout=C(1:nargout);
end