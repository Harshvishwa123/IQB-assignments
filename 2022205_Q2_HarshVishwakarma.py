import numpy as np
import pandas as pd

rx_final=[]      # List to store final alignment sequences for sequence x
ry_final=[]      # List to store final alignment sequences for sequence y

def nw(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)   # Length of sequence x
    ny = len(y)   # Length of sequence y
    
    # Initialization of the matrix.
    # Create a matrix filled with zeros of size (nx + 1) x (ny + 1)
    F = np.zeros((nx + 1, ny + 1))

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)

    # Initialize the first column with 3 (indicating vertical move)
    P[:, 0] = 3

    # Initialize the first row with 4 (indicating horizontal move)
    P[0, :] = 4

    # Matrix filling.
    current_max_score=-1     # Initialize current maximum score
    pointers=[]      # List to store pointers to optimal alignment positions

    for i in range(1, nx + 1):          # Iterate over sequence x     #COMPLETE THE CODE
        for j in range(1, ny + 1):      # Iterate over sequence y
            diagonal_score=0       # Initialize score for diagonal move
            vertical_score=0       # Initialize score for vertical move
            horizontal_score=0     # Initialize score for horizontal move
            default_score=0        # Initialize default score

            # Calculate scores for different moves based on match/mismatch/gap penalties
            if(x[i-1]==y[j-1]):       
                # If the characters at the current positions match
                diagonal_score=F[i-1][j-1]+match
            else:                     
                # If the characters at the current positions mismatch
                diagonal_score=F[i-1][j-1]+mismatch

            vertical_score=F[i-1][j]+gap         # Vertical move score
            horizontal_score=F[i][j-1]+gap       # Horizontal move score

            # Find the maximum score among all possible moves
            max_score=max(diagonal_score,vertical_score,horizontal_score,default_score)
            F[i][j]=max_score         # Update the score in the matrix

            # Update the pointers matrix based on the move that resulted in the maximum score
            if diagonal_score == max_score:
                P[i,j] = P[i,j] + 2
            if (max_score==current_max_score):
                newPositions=[i,j]
                pointers.append(newPositions)
            elif(max_score>current_max_score):
                current_max_score=max_score
                newPositions=[i,j]
                pointers.clear()
                pointers.append(newPositions)

    # Print scoring matrix using pandas DataFrame
    print("\nScoring Matrix:\n")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Trace back through the pointers to find alignments
    for pointer in pointers:
        i=pointer[0]
        j=pointer[1]
        x_character=''
        y_character=''
        rx = []      # List to store alignment sequence for sequence x
        ry = []      # List to store alignment sequence for sequence y

        # Trace back until reaching the edge of the matrix
        while P[i,j] !=0:   
            # If the move is diagonal   
            if P[i, j] in [2]:
                x_character=x[i-1]       # Add the character from sequence x
                y_character=y[j-1]       # Add the character from sequence y
                rx.append(x_character)
                ry.append(y_character)
                i=i-1
                j=j-1

        # Append the reversed alignment sequences to the final lists
        rx_final.append(rx)
        ry_final.append(ry)


# Function to calculate alignment score
def score(rx,ry):
    score=0
    for i in range(0,len(rx)):
        if(rx[i]==ry[i]):        # If characters match
            score+=2
        elif rx[i]=='-' or ry[i]=='-':        # If there is a gap
             score+=-3
        else:                  # If characters mismatch
             score+=-1
    return score
        

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
nw(seq1, seq2)            # Call the alignment function

print("\nOptimal Alignments with their scores:\n")
# Iterate over the alignments and print them with their scores
for i in range(0,len(rx_final)):
    rx=rx_final[i]
    ry=ry_final[i]
    rx = ''.join(rx)[::-1]       # Reverse the alignment sequences for printing
    ry = ''.join(ry)[::-1]
    print("Alignment",i+1)
    print('\n'.join([rx, ry]))      # Print alignment sequences
    print("score:",score(rx,ry))    # Print alignment score
    print("\n") 
