for j = 0:49
    test_js('latin hypercube', j);
    test_js('ss latinHypercube separatedLHParameters', j);
    test_js('ss latinHypercube clusteredParameters', j);
end
