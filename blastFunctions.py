"""
    This part of the BLAST heuristic algorithm defines functions
    that are going to be used in searching for the best match between
    query and a database.
"""

# This function generates k-mer's of defined length with indexes.
# The default k-mer length is equal to 3.
# Defaults:
# kmerLength = 3;
# W - kmer length.
# booleanIndexing is set to True if k-mer indexes are required.
def databaseKmerGeneration(proteinDatabaseString, kmerLength):
    databaseKmerDictionary = {}

    for i in range(0, len(proteinDatabaseString) - 2):
        # This part is used for kmer list generation
        # without indexes.
        kmer = proteinDatabaseString[i:kmerLength + i]

        # This part is used for dictionary generation
        # with indexes.
        if kmer in databaseKmerDictionary:
            databaseKmerDictionary[kmer]['index'].append(i + 1)
        else:
            databaseKmerDictionary[kmer] = {'index': [i+1]}

    return databaseKmerDictionary

# This function is used to generate high scoring sequence pairs (HSSPs).
# It is done by comparing query k-mers to the database k-mers and calculating
# the threshold value. If it is below the specified or default it is ignored.
def generatingHSSP(scoringMatrix, databaseKmerDictionary, queryFASTAString, kmerLength, thresholdValue):
    HSSPdictionary = {}

    for i in range(0, len(queryFASTAString)-2):
        queryKmer = queryFASTAString[i:kmerLength + i]

        for kmer in databaseKmerDictionary:
            T = 0
            for z in range(0, 3):
                T = T + int(scoringMatrix[kmer[z:z+1].upper()][queryKmer[z:z+1].upper()])
            if (T >= thresholdValue):
                HSSPdictionary[kmer] = {'queryIndex': i+1, 'databaseIndex': databaseKmerDictionary[kmer]['index'], 'initialHSSPValue': T}

    return HSSPdictionary


# This function generates the initial scoring matrix for the
# Smith-Waterman algorithm. Each cell of the matrix is filled with
# a dictionary, that holds the score key with value 0, and a direction
# key with a value, that is later going to be used for determining
# local alignment.
def scoreBoardInitialization(columns, rows):
    scoreBoard = []
    length = min(columns, rows)

    for i in range(length + 1):
        scoreBoard.append([{'score': 0, 'direction': 'STOP'}])
        for j in range(length):
            scoreBoard[i].append({'score': 0, 'direction': 'STOP'})
    return scoreBoard


# This function is used to return a matrix fragment. The scoring matrix
# for the Smith-Waterman algorithm is filled by 2x2 cells, not by calculating
# for example each column or row first. This means that each value for the cell
# is calculated from a 2x2 matrix which is sufficient:
# [0, 0]
# [0, x]
# x in this example is the cell that is going to be calculated.
def getMatrixFragment(scoreBoard, coordinateArray):
    matrixFragment = [[0, 0], [0, 'x']]
    matrixFragment[0][0] = scoreBoard[coordinateArray[0]-1][coordinateArray[1]-1]
    matrixFragment[0][1] = scoreBoard[coordinateArray[0]-1][coordinateArray[1]]
    matrixFragment[1][0] = scoreBoard[coordinateArray[0]][coordinateArray[1]-1]
    return matrixFragment


# This function is used to calculate the value for the specified cell.
# It is done by parsing a matrix fragment described above.
# Then a regular Smith-Waterman algorithm is run for the 2x2 matrix.
# The result is outputed as a dictionary, where score, direction for determining
# local alignment, database and query aminoacids are stored.
def smithWatermanAlgorithm(matrixFragment, scoringMatrix, gapScore, aminoAcidArray):

    vertical = matrixFragment[0][1]['score'] + gapScore
    horizontal = matrixFragment[1][0]['score'] + gapScore
    diagonal = matrixFragment[0][0]['score'] + int(scoringMatrix[aminoAcidArray[0]][aminoAcidArray[1]])

    maxScore = max(vertical, horizontal, diagonal, 0)
    if maxScore == vertical:
        return {'score': vertical, 'direction': 'vertical', 'database': '-', 'query': aminoAcidArray[1]}
    elif maxScore == horizontal:
        return {'score': horizontal, 'direction': 'horizontal', 'database': aminoAcidArray[0], 'query': '-'}
    elif maxScore == diagonal:
        return {'score': diagonal, 'direction': 'diagonal', 'database': aminoAcidArray[0], 'query': aminoAcidArray[1]}
    else:
        return {'score': 0, 'direction': 'STOP'}


# This is the function for processing every data in this program and outputing the final
# dictionary with the results.
def searchDatabaseAlgorithm(proteinDatabaseEndIndices, HSSPdictionary, proteinDatabaseDictionary, scoringMatrix, proteinDatabaseString, proteinQueryString, extensionThreshold = 14, gapScore = -10, kmerLength = 3):

    resultDictionary = {}
    dictionaryIndex = 0

    for kmer in HSSPdictionary:
        for dataIndex in HSSPdictionary[kmer]['databaseIndex']:
            databaseStartIndex = dataIndex
            databaseEndIndex = 0

            for endIndex in proteinDatabaseEndIndices:
                if endIndex == databaseStartIndex:
                    continue
                elif endIndex < databaseStartIndex:
                    continue
                else:
                    databaseEndIndex = endIndex
                    break

            if (databaseEndIndex - databaseStartIndex <= kmerLength):
                continue

            currentSequenceHeader = ''
            for header in proteinDatabaseDictionary:
                if proteinDatabaseDictionary[header]['end'] == databaseEndIndex:
                    currentSequenceHeader = header
                    break

            trueDatabase = proteinDatabaseString[dataIndex:databaseEndIndex]
            trueQuery = proteinQueryString[HSSPdictionary[kmer]['queryIndex']:]
            queryLength = len(trueQuery)
            databaseLength = len(trueDatabase)
            scoreBoard = scoreBoardInitialization(databaseLength, queryLength)

            indexRange = queryLength + 1 if databaseLength >= queryLength else databaseLength + 1

            count = 0
            maxCoordinates = [0, 0]
            lastMax = 0
            for index in range(1, indexRange):
                maxScore = 0
                aScore, bScore, cScore = {'score': 0}, {'score': 0}, {'score': 0}
                coordinateArray = [i for i in range(1, index + 1)]
                lastElement = coordinateArray.pop()
                for element in coordinateArray:
                    a = getMatrixFragment(scoreBoard, [element, lastElement])
                    b = getMatrixFragment(scoreBoard, [lastElement, element])
                    aScore = smithWatermanAlgorithm(a, scoringMatrix, gapScore, [trueDatabase[element - 1], trueQuery[lastElement - 1]])
                    bScore = smithWatermanAlgorithm(b, scoringMatrix, gapScore, [trueDatabase[lastElement - 1], trueQuery[element - 1]])
                    scoreBoard[element][lastElement] = aScore
                    scoreBoard[lastElement][element] = bScore

                    if (index == indexRange-1):
                        if (aScore['score'] > lastMax):
                            lastMax = aScore['score']
                            maxCoordinates = [element, lastElement]
                        if (bScore['score'] > lastMax):
                            lastMax = bScore['score']
                            maxCoordinates = [lastElement, element]

                c = getMatrixFragment(scoreBoard, [lastElement, lastElement])
                cScore = smithWatermanAlgorithm(c, scoringMatrix, gapScore, [trueDatabase[lastElement - 1], trueQuery[lastElement - 1]])
                scoreBoard[lastElement][lastElement] = cScore

                if (index == indexRange-1):
                    if (cScore['score'] > lastMax):
                        lastMax = cScore['score']
                        maxCoordinates = [lastElement, lastElement]

                maxScore = max(aScore['score'], bScore['score'], cScore['score'])
                count = count + 1

                if (maxScore < extensionThreshold and count >= 3):
                    break



            alignedDatabase = ''
            alignedQuery = ''
            seperator = ''
            globalMaxScore = 0
            currentCoordinates = maxCoordinates

            for reverseIndex in range(0, indexRange-1):
                currentCell = scoreBoard[currentCoordinates[0]][currentCoordinates[1]]
                direction = currentCell['direction']
                globalMaxScore = currentCell['score'] + globalMaxScore

                if (direction == 'STOP'):
                    break

                alignedDatabase = currentCell['database'] + alignedDatabase
                alignedQuery = currentCell['query'] + alignedQuery

                if (currentCell['database'] == currentCell['query']):
                    seperator = '|' + seperator
                else:
                    seperator = ' ' + seperator

                if (direction == 'diagonal'):
                    currentCoordinates = [currentCoordinates[0]-1, currentCoordinates[1]-1]
                elif (direction == 'vertical'):
                    currentCoordinates = [currentCoordinates[0]-1, currentCoordinates[1]]
                elif (direction == 'horizontal'):
                    currentCoordinates = [currentCoordinates[0], currentCoordinates[1]-1]
                else:
                    break

            resultDictionary[dictionaryIndex] = {'Header': currentSequenceHeader, 'Score': globalMaxScore, 'AlignedDatabase': alignedDatabase, 'AlignedQuery': alignedQuery, 'Seperator': seperator}
            dictionaryIndex = dictionaryIndex + 1

    return resultDictionary

# This function returns a systematic represntation
# of the results.
def displayBlastResults(sortedBlastResultsDictionary):
    alreadyPrinted = []
    for result in sortedBlastResultsDictionary:
        existBoolean = False
        object = sortedBlastResultsDictionary[result]
        if (object['Score'] == 0):
            continue
        for header in alreadyPrinted:
            if (object['Header'] == header):
                existBoolean = True
        if (existBoolean == True):
            continue
        alreadyPrinted.append(object['Header'])

        print(object['Header'] + ' with score of:', object['Score'])
        print('Database hit: ' + object['AlignedDatabase'])
        print('              ' + object['Seperator'])
        print('Query hit:    ' + object['AlignedQuery'], '\n')
