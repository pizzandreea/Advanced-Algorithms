import math
import random
import copy

f = open("output.txt", "w")

pop = int(input("Write the size of the population: "))

print("Write the function domain: ")
a = int(input("lower bound = "))
b = int(input("upped bound = "))

print("Write the Polynomial coefficients: ")
s = input().split()
poli = [int(x) for x in s]

precision = int(input("Write the precision: "))
crossoverp = float(input("Write the corssover probability: "))
mutationp = float(input("Write the mutation probability: "))
nostages = int(input("Write the number of stages: ")) + 1

# calculam lungimea cromozomului in functie de interval
lenCh = math.ceil(math.log2((b - a) * (10 ** precision)))

# generam aleator genele chromozomilor din prima populatie (pop)
# vecgtor de vectori chromosomes[pop][lenCh]
# chromosomes[0] - primul chromozom
# chromosomes[0][0] - prima gena din primul chromozom

Chromosomes = [[random.randint(0, 1) for x in range(lenCh)] for y in range(pop)]


# print(pop)

# calculul polinomului
def Fx(x, v):
    return v[0] * (x ** 2) + v[1] * x + v[2]


def Intervals(intervals, selection, stage, pop):
    # intervals = [0]
    sumi = selection[0]
    intervals.append(sumi)
    if stage == 1:
        f.write("0 " + str(sumi) + " ")
    for i in range(1, pop):
        # adaugam probabilitatea curenta
        sumi += selection[i]
        intervals.append(sumi)
        if stage == 1:
            f.write(str(sumi) + " ")
    f.write("\n\n")
    return sumi


def BinaryS(x, v, l, r):
    global last
    while l <= r:
        mid = (l + r) // 2
        if v[mid] <= x:
            last = mid
            l = mid + 1
        elif v[mid] > x:
            r = mid - 1
    return last + 1


for stage in range(1, nostages):
    if stage == 1:
        f.write("First Population\n")

    idFittest = 0

    maxFittness = float('-inf')
    # all x's
    Xs = []

    # sum of all F
    sumF = 0

    # f.write(" %s " % str(Chromosomes[1]))
    # pentru fiecare cromozom din populatia curenta
    for i in range(pop):
        ############### calcularea f(X) si cautam cel mai mare valoare a functiei ##############################

        stringCh = ''.join([str(x) for x in Chromosomes[i]])
        x = int(stringCh, 2)  # transform in baza 10
        # codificarea unui cromozom - intrapolare pe D - formula
        intrapolateX = ((b - a) / (2 ** lenCh - 1)) * x + a
        # precision zecimals
        roundIntrapolateX = round(intrapolateX, precision)

        # calculez F(x)
        fX = Fx(roundIntrapolateX, poli)

        # memorez toti x cromozomilor
        Xs.append(intrapolateX)
        sumF += Fx(intrapolateX, poli)

        if fX > maxFittness:
            maxFittness = fX
            idFittest = i

        if stage == 1:
            f.write("    %s: %s x= %s f= %s\n" % (str(i + 1), stringCh, str(roundIntrapolateX), str(fX)))

    fittest = Chromosomes[idFittest].copy()

    ##################### Calculam probabilitatile de selectie ######################
    if stage == 1:
        f.write("\nSelection Probabilities\n")
    selProb = []
    for i in range(pop):
        prob = Fx(Xs[i], poli) / sumF
        selProb.append(prob)
        if stage == 1:
            f.write("cromozom %s probability %s\n" % (str(i + 1), str(prob)))

    ######################## Aflam intervalele de selectie ############################
    if stage == 1:
        f.write("\nSelection probabilities Intervals\n")
    #     intervalele vor fi suma probabilitatilor
    selProbIntervals = [0]

    sumI = Intervals(selProbIntervals, selProb, stage, pop)
    # f.write(str(selProbIntervals[1]))
    ####################### Randomly selecting chromosomes #############################

    selected = [0 for x in range(pop)]

    for i in range(pop):
        r = random.random()
        randCh = BinaryS(r, selProbIntervals, 0, pop) - 1
        if stage == 1:
            f.write("u = %s selected chromosome: %s\n" % (r, randCh))
        selected[i] = randCh

    if stage == 1:
        f.write("\nAfter Selection\n")
    selectedCh = []
    for i in range(pop):
        if stage == 1:
            bit = ''.join([str(x) for x in Chromosomes[selected[i]]])
            fX = Fx(Xs[selected[i]], poli)
            f.write("%s: %s x= %s f=%s\n" % (str(i+1), bit, str(round(Xs[selected[i]], precision)), fX))
        selectedCh.append(Chromosomes[selected[i]])

    # keeping the newly selected chromosomes for the next stage
    Chromosomes = copy.deepcopy(selectedCh)

    ########################### Crossover ################################

    # lista de indici
    crossovers = []
    if stage == 1:
        f.write("\nCrossover Probability %s\n" % (crossoverp))

    for i in range(pop):
        r = random.random()
        bit = ''.join([str(x) for x in Chromosomes[selected[i]]])

        if stage == 1:
            f.write("%s: %s u=%s" % (str(i + 1), bit, str(r)))
            if r < crossoverp:
                f.write(" < %s participates in crossover\n" % (crossoverp))
                crossovers.append(i)
            else:
                f.write("\n")

    if stage == 1:
        f.write("\n")

    while len(crossovers) > 1:
        i = len(crossovers)-1
        j = i - 1

        if stage == 1:
            f.write("\nCrossover between chromosome %s and chromosome %s: " % (str(crossovers[i]+1), str(crossovers[j]+1)))
        # the point where we slice the chromosomes
        slice = random.randrange(lenCh)

        if stage == 1:
            bit1 = ''.join([str(x) for x in Chromosomes[crossovers[i]]])
            bit2 = ''.join([str(x) for x in Chromosomes[crossovers[j]]])
            f.write("\n%s and %s slicing point: %s" % (bit1, bit2, str(slice)))
            f.write("\n")

        # we switch the first 'slice' digits of the chromosome a
        first = Chromosomes[crossovers[i]][:slice+1].copy()
        Chromosomes[crossovers[i]][:slice+1] = Chromosomes[crossovers[j]][:slice+1].copy()
        Chromosomes[crossovers[j]][:slice+1] = first.copy()

        bit1 = ''.join([str(x) for x in Chromosomes[crossovers[i]]])
        bit2 = ''.join([str(x) for x in Chromosomes[crossovers[j]]])

        if stage == 1:
            f.write("Result: + %s and %s \n" % (bit1, bit2))

        aux = crossovers[:len(crossovers)-1]
        crossovers = aux.copy()

    if stage == 1:
        f.write("\nAfter Crossover\n")

    for i in range(pop):
        bit = ''.join([str(x) for x in Chromosomes[i]])
        x = int(bit,2)
        intrapolateX = ((b - a) / (2 ** lenCh - 1)) * x + a
        Xs[i] = intrapolateX
        roundIntrapolateX = round(intrapolateX, precision)
        fX = Fx(roundIntrapolateX, poli)
        if stage == 1:
            f.write("%s : %s x= %s f= %s\n"  % (str(i+1), bit, roundIntrapolateX, fX))


    if stage == 1:
        f.write("\nThe mutation probability for every gene is %s\n" % (str(mutationp)))
        f.write("The following chromosomes have mutated: \n")
    for i in range(pop):
        r = random.random()
        if r < mutationp:
            gene = random.randrange(lenCh)
            # the gene that mutates switches from 1 to 0
            Chromosomes[i][gene] = 1-Chromosomes[i][gene]
            if stage == 1:
                f.write(str(i+1) +"\n")

    if stage == 1:
        f.write("\n After mutation:\n")
#     se repeta afisarile de mai sus
