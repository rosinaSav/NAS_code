import time
import os
import re
import numpy as np

def get_datasets():

    '''
    Extact the ESEs from the data files
    '''

    ese_motifs = {}

    for filePath in os.listdir('data/'):
        with open('data/{}'.format(filePath), 'rU') as openFile:
            lines = openFile.readlines()

            dataset = lines[0].strip('\n')
            ese_motifs[dataset] = []

            for line in lines[1:]:
                ese_motifs[dataset].append(line.strip('\n'))

    return(ese_motifs)



def getNtPairs(eses):

    nt_pairs = []
    for ese in eses:
        nt_pairs.extend(re.findall(r'.{2}', ese))

    return(nt_pairs)


def generateRandomMotifs(eseSet, nt_pairs):

    randomMotifs = []
    while len(randomMotifs) < len(eseSet):
        randomMotif = ''.join(np.random.choice(nt_pairs, 3))
        if randomMotif not in eseSet and randomMotif not in randomMotifs:
            randomMotifs.append(randomMotif)

    return(randomMotifs)

def countStops(list):

    count = 0
    for item in list:
        count += len(re.findall(r'TGA|TAA|TAG', item))

    return(count)

def getLongestMononucleotideRun(set):

    longestItem = max(len(i) for i in set)
    print(longestItem)

def runData(data):

    dataCount = 0

    for dataset in data:
        dataCount += 1
        eseSet = data[dataset]

        if dataCount <= 1:
            real_count = countStops(eseSet)
            getLongestMononucleotideRun(eseSet)

            randomisedCounts = []

            nt_pairs = getNtPairs(eseSet)

            for i in range(20):
                randomMotifs = generateRandomMotifs(eseSet, nt_pairs)
                randomisedCounts.append(countStops(randomMotifs))

            print(real_count)
            print(randomisedCounts)

def run():

    data = get_datasets()
    runData(data)




def main():

    t0_main = time.time()
    run()
    t1_main = time.time()
    print ('\n%s\nTime: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
