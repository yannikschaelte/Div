clear all;

for j = 0:49
    test_rme('latin hypercube', j);
    test_rme('ss latinHypercube separatedLHParameters', j);
    test_rme('ss latinHypercube clusteredParameters', j);
end
