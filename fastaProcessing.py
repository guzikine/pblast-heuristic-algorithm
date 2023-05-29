'''
    This part of the BLAST heuristic algorithm
    extracts seq data from fasta file as well as protein
    coding positions in the database string.
'''

import regex

# booleanDictionary is set to TRUE if dictionary output is required.
# idDictionary returns a dictionary variable that can be accesses as follows
# as an example
#
# dictionaryValue = idDictionary["XP_009096760.1 albumin [Serinus canaria]"]
# -> {'start': 609, 'end': 1224}
#
# dictionaryValueEnd = idDictionary["XP_009096760.1 albumin [Serinus canaria]"]["end"]
# -> 1224
#
# dictionaryValueStart = idDictionary["XP_009096760.1 albumin [Serinus canaria]"]["start"]
# -> 609
def readFastaFile(fastaFile, booleanDatabase = False):
    fileLines = fastaFile.split("\n")
    headerPattern = '(^>(.+))'
    fastaString = ''
    idDictionary = {}
    lengthIndex = 0
    booleanFirstSequence = True
    previousHeader = ''
    databaseEndIndices = []

    for line in fileLines:
        headerMatch = regex.match(headerPattern, line)
        if headerMatch and booleanFirstSequence:
            booleanFirstSequence = False
            previousHeader = headerMatch.group(2)
            idDictionary[headerMatch.group(2)] = {'start': lengthIndex}

        elif headerMatch:
            idDictionary[previousHeader]['end'] = lengthIndex
            idDictionary[headerMatch.group(2)] = {'start': lengthIndex}
            databaseEndIndices.append(lengthIndex)
            previousHeader = headerMatch.group(2)

        elif line == fileLines[-1]:
            idDictionary[previousHeader]['end'] = lengthIndex
            databaseEndIndices.append(lengthIndex)

        else:
            fastaString = fastaString + line
            lengthIndex = lengthIndex + len(line)

    if booleanDatabase:
        return fastaString, idDictionary, databaseEndIndices
    else:
        return fastaString