"""
    This part of the BLAST heuristic algorithm will be
    used to process scoring matrices and prepare them for
    use in the algorithm.

    BLOSUM matrices:
    X - any amino acid.
    B - aspartic acid or asparagine.
    Z - glutamic acid or glutamine.
"""

import pandas as pd
import regex

# This function reads the scoring matrix and places it into a dataframe.
def readScoringMatrixFile(scoringMatrixFile):
    fileLines = scoringMatrixFile.split("\n")
    aminoAcidPattern = '((^[A-Z]{1}).+)'
    dataframeColumnsRows = []
    matrixLines = []

    for line in fileLines:
        aminoAcidMatch = regex.match(aminoAcidPattern, line)
        if aminoAcidMatch:
            matrixLines.append(aminoAcidMatch.group(1))
            dataframeColumnsRows.append(aminoAcidMatch.group(2))

    matrixDataFrame = pd.DataFrame(index=dataframeColumnsRows, columns=dataframeColumnsRows)

    scorePattern = '(^[A-Z]{1}\s*)(((.?\d)\s*)+)'
    for i in range(0, len(matrixLines)):
        scoreMatch = regex.match(scorePattern, matrixLines[i])
        for j in range(0, len(matrixLines)):
            matrixDataFrame.iloc[i ,j] = scoreMatch.captures(4)[j]

    return matrixDataFrame