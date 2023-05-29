"""
    This is the main part of the BLAST heuristic algorithm.
    This is where all the functions are executed and files are imported.
"""

import blastFunctions
import fastaProcessing
import matrixProcessing

# Defining file and library variables for passing in
# read files.
scoringMatrixFile = ''
queryFasta = ''
databaseFasta = ''

try:
    # Opening matrix file for reading.
    fMatrix = open("matrices/BLOSUM62.matrix", mode='r', encoding='utf-8')

    # Opening query FASTA file for reading.
    fQuery = open("query.fasta")

    # Opening FASTA protein databse file for reading.
    fDatabase = open("proteinDatabase.fasta")

    # Reading each file and storing it in a variable.
    scoringMatrixFile = fMatrix.read()
    queryFasta = fQuery.read()
    databaseFasta = fDatabase.read()
except Exception as e:
    print(e)
    print("Something went wrong when opening the file to read.")
finally:
    # Closing each file after reading was done.
    fMatrix.close()
    fQuery.close()
    fDatabase.close()

def runpBlast(kmerLength = 3, thresholdValue = 13, gapScore = -10):
    extensionThreshold = thresholdValue + 1

    # Passing matrix file to the matrixProcessing function that is described in
    # matrixProcessing.py. This function creates a pandas module dataFrame.
    scoringMatrix = matrixProcessing.readScoringMatrixFile(scoringMatrixFile)

    # Passing FASTA database file to the fastaProcessing function that is described in
    # fastaProcessing.py. This function converts FASTA file into a single FASTA string
    # that is going to be used as a database. The second parameter can be set to TRUE
    # to also return a dictionary that defines where each protein molecule sequence
    # starts and ends.
    proteinDatabaseString, proteinDatabaseDictionary, proteinDatabaseEndIndices = fastaProcessing.readFastaFile(databaseFasta, True)

    # Passing query FASTA file to the fastaProcessing function that is described in
    # fastaProcessing.py. This funciton converts FASTA file into a single FASTA string
    # that is going to be used as a database.
    proteinQueryString = fastaProcessing.readFastaFile(queryFasta)

    # Generating FASTA database sequence k-mers with their indexes.
    databaseKmerDictionary = blastFunctions.databaseKmerGeneration(proteinDatabaseString, kmerLength)

    # Generating high scoring sequence pairs between the database and the query.
    HSSPdictionary = blastFunctions.generatingHSSP(scoringMatrix, databaseKmerDictionary, proteinQueryString, kmerLength, thresholdValue)

    # This function returns a dictionary with the results after running the actual pBLAST algorithm.
    blastResultsDictionary = blastFunctions.searchDatabaseAlgorithm(proteinDatabaseEndIndices, HSSPdictionary, proteinDatabaseDictionary, scoringMatrix, proteinDatabaseString, proteinQueryString, extensionThreshold, gapScore, kmerLength)

    # This function sorts the resultDictionary in ascending manner by the ['Score'] value.
    sortedBlastResultsDictionary = dict(sorted(blastResultsDictionary.items(), key=lambda item: item[1]['Score'], reverse=True))

    blastFunctions.displayBlastResults(sortedBlastResultsDictionary)