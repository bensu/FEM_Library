function index = index_range(step,counter,varargin)
% index = index_range(step,counter)
% index = index_range(step,counter,index2)
    aux = ((step*(counter-1)+1):(step*counter))';
    if isempty(varargin)
        index = aux;
    else index = aux(varargin{1});
    end
end