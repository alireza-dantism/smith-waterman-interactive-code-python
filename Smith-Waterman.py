# -*- coding: utf-8 -*-
# Smith Waterman algorithm as an interactive code using Python - Using 1 Matrix
# By Alireza Dantism
# Created on Thursday, June 30, 2022
# Edited on July 2022

# Access to the arguments passed to the script
import sys

# Function for reading files and get sequences
def read_sequence_file(filename):
    filereader = open(filename)
    line = filereader.readline()
    line = line.replace("\n", "")
    sequence = ""
    while True:
        line = filereader.readline()
        if not line:
            break
        sequence += line.replace("\n", "")
    return sequence

# Function for showing help
def helpApp():
    print("\n \n \n Hi Dear user, welcome to the program  ")
    print(" Here is some user guid for you:      \n")
    print(" ----------------------")
    print(" ****** python smith.py file1.fasta file2.fasta -A S/N -M 1 -m -2 -g -3   ****** ")
    print(" ----------------------")
    print("\n You can pass below parameter to code ")
    print("\n python smith.py file1.fasta file2.fasta => for passing Fasta files   ")
    print("\n python smith.py -A S/N => for passing algorithm Smith-Waterman or Needleman-Wunsch    ")
    print("\n python smith.py -M 1   => for passing match score                  ")
    print("\n python smith.py -m -1  => for passing missmatch score              ")
    print("\n python smith.py -g -2  => for passing gap score                    ")
    print("\n default values are =>  \"TCAGTTGCC\" \"AGGTTG\" match = 1, missmatch = -1 , gap score = -2  ******  ")
    input()


# Global Variables with default values
ALGORITHM_TYPE = "S"
MATCH_VALUE = 1
MISS_MATCH_VALUE = -1
GAP_VALUE = -2
parameter_1 = ""
sequence1 = "TCAGTTGCC"
sequence2 = "AGGTTG"
seq_check = 0

# Checking arguments passed to script and change default values
for x in range(1, len(sys.argv)):
    # -A => Algorithm Type
    if sys.argv[x] == "-A":
        ALGORITHM_TYPE = sys.argv[x+1]
    # -M => Match Score
    elif sys.argv[x] == "-M":
        MATCH_VALUE = int(sys.argv[x+1])
    # -m => Missmatch Score
    elif sys.argv[x] == "-m":
        MISS_MATCH_VALUE = int(sys.argv[x+1])
    # -g => Gap Penalty
    elif sys.argv[x] == "-g":
        GAP_VALUE = int(sys.argv[x+1])
    # Getting Fasta Files Name
    elif (".fasta" in sys.argv[x] or ".fa" in sys.argv[x]) and (seq_check == 0):
        sequence1 = read_sequence_file(sys.argv[x])
        seq_check = 1
    elif (".fasta" in sys.argv[x] or ".fa" in sys.argv[x]) and seq_check:
        sequence2 = read_sequence_file(sys.argv[x])
    # -h or -help
    elif sys.argv[x] == "-h" or sys.argv[x] == "-help":
        helpApp()


if (len(sys.argv) == 1):
    helpApp()


# 2 dimensional array for calculating Smith–Waterman scores
Matrix       = {}
Matrix[0, 0] = ""
Matrix[0, 1] = "_"

# At the end of algorithm - to store maximum value and it's coordinate we use a variable as below.
L_C_V = 0 # Last Char Value
L_C_A = 0 # Last Char Adderss

# I use this variable for putting longer string length at the top of Matrix.
flag_seq1 = sequence1

print("\n ****** Hi Dear user, welcome to the program ******  \n")

# Check which sequence has a more amino acid.
if len(sequence1) > len(sequence2):
  flag    = len(sequence1)
  Max_STR = len(sequence1)
  Min_STR = len(sequence2)
else:
  flag    = len(sequence2)
  Max_STR = len(sequence2)
  Min_STR = len(sequence1)
  sequence1 = sequence2
  sequence2 = flag_seq1

# Initialize - first row and column of Matrix score.
# In first row and column of value and coordinate matrix we show characters.

for x in range(len(sequence1)):
  Matrix[0, x+2] = sequence1[x]

Matrix[0, 0] = "0"
Matrix[1, 0] = "_"
for z in range(len(sequence2)):
  Matrix[z+2, 0] = sequence2[z]

# initialize - second row and column of Matrix score
# in second row & column of value and coordinate matrix we show zero value -but for Needleman-Wunsch we should calculate

for x in range(len(sequence1) + 2):
    if ALGORITHM_TYPE == "S":
        # Smith-Waterman
        Matrix[1, x] = 0
    else:
        # Needleman-Wunsch
        if x == 0 or x == 1:
            Matrix[1, x] = 0
        else:
            Matrix[1, x] = Matrix[1, x-1] + GAP_VALUE

for x in range(len(sequence2) + 2):
    if ALGORITHM_TYPE == "S":
        # Smith-Waterman
        Matrix[x, 1] = 0
    else:
        # Needleman-Wunsch
        if x == 0 or x == 1:
            Matrix[x, 1] = 0
        else:
            Matrix[x, 1] = Matrix[x-1, 1] + GAP_VALUE


# here Smith-Waterman algorithm start and calculate each cell value.
for i in range(len(sequence1)):
    for j in range(len(sequence2)):
        
        # MATCH
        if sequence1[i] == sequence2[j]:
            M = Matrix[j+1, i+1] + MATCH_VALUE
            m = -2000000
        # MIS MATCH
        else:
            m = Matrix[j+1, i+1] + MISS_MATCH_VALUE
            M = -2000000
               
        # Gap Penalty
        gTop  = Matrix[j+2, i+1] + GAP_VALUE
        gLeft = Matrix[j+1, i+2] + GAP_VALUE
        
        # Get max value from scores.
        score_list = [M, m, gTop, gLeft]
        max_value = max(score_list)
        max_index = score_list.index(max_value)       
            
        # Check if max value refers to negative value then replace it with zero value.
        if ALGORITHM_TYPE == "S":
            if max_value > 0:
                Matrix[j+2, i+2] = max_value
            else:
                Matrix[j+2, i+2] = 0
        else:
            Matrix[j+2, i+2] = max_value


print("\n")
#print("--------------- Here Is Smith-Waterman Matrix ---------------")
#print("\n")

# We use temporary variable for check the length of sequences
if len(sequence1) < len(sequence2):
    a = Max_STR
    b = Min_STR
else:
    a = Min_STR
    b = Max_STR

# Here we plus 2 because first and second of rows and column is amino acid characters and zero value
for i in range(a + 2):
    for j in range(b + 2):
        print(Matrix[i, j], end=" \t")
        if i >= 2 and j >= 2:
            if L_C_V < Matrix[i, j]:
                L_C_V = Matrix[i, j]
                L_C_A = str(i) + ',' + str(j)
    print("\n")
print("\n")

#print("** Max value in matrix is ", L_C_V, " in coorinates of ", "Matrix[" + L_C_A + "] \n")

# Generate aligned sequence using these variables.
firstSeq = ""
secondSeq = ""
flagR = 0
i = 1

while i == 1:

    row = int(L_C_A.split(',')[0])
    col = int(L_C_A.split(',')[1])

    if (row == 1) and (col == 1):
        i = 0
        break

    # Get amino acid character
    AminoAcidTop  = str(Matrix[row - row, col])
    AminoAcidLeft = str(Matrix[row, col - col])

    # Get cells values
    ThisValue = Matrix[row, col]
    LeftValue = Matrix[row, col - 1]
    TopValue  = Matrix[row - 1, col]
    SkewValue = Matrix[row - 1, col - 1]

    if ThisValue == 0 and ALGORITHM_TYPE == "S":
        i = 0
        break

    if not (AminoAcidTop == "0" or AminoAcidLeft == "0"):
        # Sort Values => For Checking Index
        x_s = sorted([LeftValue, TopValue, SkewValue])
    else:
        # NeedleMan
        if AminoAcidLeft == "0":
            firstSeq += "*"
            secondSeq += AminoAcidTop

            if row - 1 == 0:
                L_C_A = str(row) + ',' + str(col - 1)
                continue

        elif AminoAcidTop == "0":
            firstSeq += AminoAcidLeft
            secondSeq += "-"

            if col - 1 == 0:
                L_C_A = str(row - 1) + ',' + str(col)
                continue

        if LeftValue == 0 and TopValue == 0 and ThisValue == 0:
            i = 0
            break

    if (x_s[2] - x_s[1]) == 1:
        if (x_s[1] - x_s[0]) == 1:
            if ALGORITHM_TYPE == "N":
                flagR = 1
        else:
            flagR = 0
    else:
        flagR = 0

    max_value = max(LeftValue, TopValue, SkewValue)


    if (flagR == 1) or ((ALGORITHM_TYPE == "S") and (TopValue == LeftValue == SkewValue)):
        firstSeq += AminoAcidLeft
        secondSeq += AminoAcidTop
        L_C_A = str(row - 1) + ',' + str(col - 1)
    else:
        if AminoAcidLeft == AminoAcidTop:
            firstSeq += AminoAcidLeft
            secondSeq += AminoAcidTop
            L_C_A = str(row - 1) + ',' + str(col - 1)
        else:
            if max_value == LeftValue:
                firstSeq += "-"
                secondSeq += AminoAcidTop
                L_C_A = str(row) + ',' + str(col - 1)
            elif max_value == TopValue:
                firstSeq += AminoAcidLeft
                secondSeq += "-"
                L_C_A = str(row - 1) + ',' + str(col)
            else:
                firstSeq += AminoAcidLeft
                secondSeq += AminoAcidTop
                L_C_A = str(row - 1) + ',' + str(col - 1)

    if ALGORITHM_TYPE == "S":
        if max_value == 0:
            i = 0


# reverse characters and print
print(secondSeq[::-1])
print(firstSeq[::-1])

input("\n -- Press Enter key to exit.")