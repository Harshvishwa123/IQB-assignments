import numpy as np
import pandas as pd

# Lists to store final alignments
rx_final = []
ry_final = []

# Function for traceback to find all optimal alignments
def traceback(i, j, P, rx, ry, x, y):
    # i (int): Current row index.
    # j (int): Current column index.
    # Base case: if both pointers reach the beginning of sequences
    if (i <= 0 and j <= 0):
        newx = list.copy(rx)
        newy = list.copy(ry)
        # Append the current alignment to the final lists
        rx_final.append(newx)
        ry_final.append(newy)
        return
    
    # Check the pointers and move accordingly based on the values in matrix
    if P[i, j] in [2, 5, 6, 9]:           # COMPLETE THE CODE
        # Move diagonally (match/mismatch)
        x_character = x[i-1]
        y_character = y[j-1]
        rx.append(x_character)
        ry.append(y_character)
        traceback(i-1, j-1, P, rx, ry, x, y)
        rx.pop()
        ry.pop()

    if P[i, j] in [3, 5, 7, 9]:           # COMPLETE THE CODE
        # Move vertically (gap in sequence x)
        x_character = x[i-1]
        y_character = '-'
        rx.append(x_character)
        ry.append(y_character)
        traceback(i-1, j, P, rx, ry, x, y)
        rx.pop()
        ry.pop()

    if P[i, j] in [4, 6, 7, 9]:           # COMPLETE THE CODE
        # Move horizontally (gap in sequence y)
        x_character = '-'
        y_character = y[j-1]
        rx.append(x_character)
        ry.append(y_character)
        traceback(i, j-1, P, rx, ry, x, y)
        rx.pop()
        ry.pop()

# Function to calculate score of an alignment
def score(rx, ry,match,mismatch,gap):
    score = 0
    for i in range(0, len(rx)):
        if (rx[i] == ry[i]):        # If characters match
            score += match         
        elif rx[i] == '-' or ry[i] == '-':    # If there is a gap
            score += gap                # Gap penalty should be added here
        else:         # If characters mismatch
            score += mismatch
    return score

# Function for Needleman-Wunsch algorithm
def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)

    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)    # Corrected: Gap penalty initialization.
    F[0, :] = np.linspace(0, gap * ny, ny + 1)    # Corrected: Gap penalty initialization.

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3         # Set pointers for vertical movement in the first column
    P[0, :] = 4         # Set pointers for horizontal movement in the first row

    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):  # COMPLETE THE CODE
            diagonal_score=0       # Initialize score for diagonal move
            vertical_score=0       # Initialize score for vertical move
            horizontal_score=0     # Initialize score for horizontal move

            if (x[i-1] == y[j-1]):
                diagonal_score = F[i-1][j-1]+match
            else:
                diagonal_score = F[i-1][j-1]+mismatch

            vertical_score = F[i-1][j]+gap
            horizontal_score = F[i][j-1]+gap
            max_score = max(diagonal_score, vertical_score, horizontal_score)

            F[i][j] = max_score
            if diagonal_score == max_score:
                P[i, j] = P[i, j] + 2
            if vertical_score == max_score:
                P[i, j] = P[i, j] + 3
            if horizontal_score == max_score:
                P[i, j] = P[i, j] + 4

    # Print scoring matrix using pandas DataFrame
    print("\nScoring Matrix:\n")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Trace through an optimal alignment.
    traceback(nx, ny, P, [], [], x, y)

# DNA sequences to be used:
seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
# Calling Needleman-Wunsch algorithm function
nw(seq1, seq2)

# Printing optimal alignments and their scores
print("\nOptimal Alignments with their scores:\n")
for i in range(0, len(rx_final)):
    rx = rx_final[i]
    ry = ry_final[i]
    rx = ''.join(rx)[::-1]       # Reverse the alignment sequences for printing
    ry = ''.join(ry)[::-1]
    print("Alignment",i+1)
    print('\n'.join([rx, ry]))    # Print alignment sequences
    print("Score: ", score(rx, ry,match=2,mismatch=-3,gap=-1))
    print("\n")
