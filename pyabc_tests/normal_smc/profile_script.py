import pstats

p = pstats.Stats("tmp.dat")
p.sort_stats('cumulative')
p.print_stats()
