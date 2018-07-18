clear all;

for j = 0:49
    disp(['Random seed ' num2str(j)]);
    
    test_js('latin hypercube', j);
    test_js('ss latinHypercube separatedLHParameters', j);
    test_js('ss latinHypercube clusteredParameters', j);
end
