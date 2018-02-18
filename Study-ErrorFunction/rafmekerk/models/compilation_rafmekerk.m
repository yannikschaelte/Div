[exdir,~,~]=fileparts(which('rafmekerk_normal_standard_syms.m'));

amiwrap('rafmekerk_normal_standard','rafmekerk_normal_standard_syms',exdir)
amiwrap('rafmekerk_laplace_standard','rafmekerk_laplace_standard_syms',exdir)
amiwrap('rafmekerk_generalizednormal_standard','rafmekerk_generalized_normal_standard_syms',exdir)