
out = open('data/tennessen_popsize_fine_X.txt', 'w')

with open('data/tennessen_popsize_fine.txt') as f:
    for line in f:
        t, p = line.split()
        p = int(int(p) / 2)
        print(t, p, sep='\t', file=out)


