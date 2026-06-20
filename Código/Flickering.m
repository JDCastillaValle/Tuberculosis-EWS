
% Dado un arreglo (array) y un valor (value), la función dice 
% cuantas veces el "valor" queda entre dos valores consecutivos del arreglo 
% recorriendolo de forma ascendente. Y regresa una lista con los índices 
% de los pares de valores del arreglo que contienen al valor.

function [flicks, index] = Flickering(array, value)

    a = array(1:end-1);
    b = array(2:end);

    between = (value >= min(a, b)) & (value <= max(a, b));

    index = find(between);
    flicks = length(index);
end
