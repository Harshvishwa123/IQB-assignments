# Bioinformatics and Structural Biology Assignments

This repository contains Python scripts submitted as part of coursework for INTRODUCTION TO QUANTITATIVE BIOLOGY. Each script tackles a distinct problem, including sequence alignment and protein secondary structure prediction.

---

## 📁 Files

### 1. `2022205_Q1_HarshVishwakarma.py`

**Topic:** Global Sequence Alignment using Needleman-Wunsch Algorithm

- Implements the classic Needleman-Wunsch algorithm to find **all optimal alignments** between two DNA sequences.
- Includes:
  - Dynamic programming for matrix filling.
  - Traceback function to explore multiple optimal paths.
  - Scoring system with configurable match, mismatch, and gap penalties.
- Outputs all alignments and their scores.

**Run:**
python 2022205_Q1_HarshVishwakarma.py

### 2. `2022205_Q2_HarshVishwakarma.py`

**Topic:** Local Sequence Alignment (Smith-Waterman variant)

- Identifies optimal local alignments between DNA sequences.
- Records highest-scoring subsequences based on:
  - Match, mismatch, and gap scores.
  - Traceback from local maxima in the matrix.
- Outputs alignments and scores.

**Run:**
bash
python 2022205_Q2_HarshVishwakarma.py

## 3. `IQB_assignment_2.py`

**Topic:** Secondary Structure Prediction using Chou-Fasman Algorithm

- Analyzes a given protein sequence to predict:
  - Alpha helices (H)
  - Beta strands (S)
- Uses:
  - Amino acid propensities from Chou-Fasman parameters.
  - Sliding window approach for nucleation detection.
  - Extension logic for elongating secondary structures.
  - Conflict resolution when both structures are predicted in overlapping regions.
- Provides a final annotated sequence indicating secondary structural regions.

**Run:**
bash
python IQB_assignment_2.py


