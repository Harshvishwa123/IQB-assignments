# Given protein sequence
protein_sequence=(  "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAP"
                    "STVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPA"
                    "KSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKT"
                    "HMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTAS"
                    "ELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRH"
                    "LRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDER"
                    "PYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCD"
                    "VRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRV"
                    "HSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVT"
                    "QRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLE"
                    "PSQDL"   )

# Coversion of One letter Code to Three letter code of different Amino Acids
One_letterCODE_to_Three_letterCODE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp',
    'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly',
    'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
    'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

# Chou and Fasman parameters to be used for the prediction
alpha = {'Glu': 1.53, 'Ala': 1.45, 'Leu': 1.34, 'His': 1.24, 'Met': 1.20, 
         'Gln': 1.17, 'Trp': 1.14, 'Val': 1.14, 'Phe': 1.12, 'Lys': 1.07, 
         'Ile': 1.00, 'Asp': 0.98, 'Thr': 0.82, 'Ser': 0.79, 'Arg': 0.79,
         'Cys': 0.77, 'Asn': 0.73, 'Tyr': 0.61, 'Pro': 0.59, 'Gly': 0.53 }

beta={'Met': 1.67,'Val': 1.65,'Ile': 1.60,'Cys': 1.30,'Tyr': 1.29,
      'Phe': 1.28,'Gln': 1.23,'Leu': 1.22,'Thr': 1.20,'Trp': 1.19,
      'Ala': 0.97,'Arg': 0.90,'Gly': 0.81,'Asp': 0.80,'Lys': 0.74,
      'Ser': 0.72,'His': 0.71,'Asn': 0.65,'Pro': 0.62,'Glu': 0.26 }

# The length of the given protein sequence
len_of_sequence=len(protein_sequence)

# Function notation that generates a string representing the secondary structure (helix or beta-strand) based on the given print data (0 or 1) for a specified range of the sequence.
def Generate_Secondary_Structure(Helix_or_BetaStrand, start, end, print_data):
    answer = ""
    char = "H" if Helix_or_BetaStrand == 1 else "S"
    i = start
    while i < start + end:
        if print_data[i] == 0:
            answer += " "
        else:
            answer += char
        i += 1
    return answer

# This function scans the entire sequence with a specified window size to identify potential nucleation sites for helices or beta-strands.
def Finding_Nucleation_Sites(window_size, Helix_or_BetaStrand, window_score):
    pointers = []
    for current_index in range(len_of_sequence):
        start_index = current_index
        end_index = current_index + (window_size - 1)
        window_sequence = protein_sequence[start_index:end_index + 1]
        if len(window_sequence) != window_size:
            break
        counter = sum(1 for char in window_sequence if One_letterCODE_to_Three_letterCODE[char] in alpha if alpha[One_letterCODE_to_Three_letterCODE[char]] >= 1)
        if counter >= window_score:
            pointers.append([start_index, end_index])
    # Return a list of pointers indicating the start and end indices of the identified nucleation sites.
    return pointers
    
# This function checks whether a given sequence of amino acids can form a valid secondary structure extension.
def Extension_for_helix_OR_strand(right_start,right_end,helix_or_BetaStrand):
    subsequence = protein_sequence[right_start:right_end + 1]
    score = sum(alpha[One_letterCODE_to_Three_letterCODE[char]] if helix_or_BetaStrand == 1 else beta[One_letterCODE_to_Three_letterCODE[char]] for char in subsequence)
    # Check the score is less than 4 or that the average P(H)<1).
    if(round(score,2)<4):
        return False
    else:
        return True

# This Function extends the identified nucleation sites for helices or beta-strands.
def extendWindow(pointers,length_of_extended_window,helix_or_betaStrand):
    final_extended_window=[]
    # It extends the identified nucleation sites to the left and right to find the maximum possible extension .
    for pointer in pointers:
        start_index=pointer[0]
        end_index=pointer[1]
        # extend this right
        old_start=start_index
        old_end=end_index
        right_start=old_end-(length_of_extended_window-1)
        right_end=old_end
        left_start=old_start
        left_end=old_start+(length_of_extended_window-1)
        while True:
            right_start=right_start+1
            right_end=right_end+1
            # check for valid string
            if(right_end>=len_of_sequence):
                old_end=right_end-1
                break
            result = Extension_for_helix_OR_strand(right_start,right_end,helix_or_betaStrand)
            if(result==True):
                continue
            else:
                old_end=right_end-1
                break

        while True:
            left_start=left_start-1
            left_end=left_end-1
            # check for valid string
            if(left_start<0):
                old_start=left_start+1
                break
            result = Extension_for_helix_OR_strand(left_start,left_end,helix_or_betaStrand)
            if(result==True):
                continue
            else:
                old_start=left_start+1
                old_end=right_end-1
                break

        final_extended_window.append([old_start,old_end])
    return final_extended_window


# Function to merge the pointers obtained from extending the nucleation sites and generate a list of 0s and 1s.
def Generate_0_1_List(final_pointers):
    print_data = []
    previous_zero_position = 0
    previous_end = -1
    for data in final_pointers:
        start, end = data
        if previous_end >= start:
            start = previous_end + 1
        if start > end:
            continue
        zero_counts = start - previous_zero_position
        one_counts = end - start + 1
        previous_zero_position = end + 1
        previous_end = end
        print_data.extend([0] * zero_counts)  # 1 represents the presence of a secondary structure element 
        print_data.extend([1] * one_counts)   # 0 represents no secondary structure

    remaining_zero = len_of_sequence - previous_zero_position
    print_data.extend([0] * remaining_zero)

    return print_data

# It prints the protein sequence along with the corresponding secondary structure elements.
def Display_The_Sequence_Regions(print_data, helix_or_betaStrand, line_length):
    for position in range(0, len_of_sequence, line_length):
        start = position
        end = position + line_length
        protein = protein_sequence[start:end]
        print()
        print(protein)
        helix_or_sheet_notation = Generate_Secondary_Structure(helix_or_betaStrand, start, len(protein), print_data)
        print(helix_or_sheet_notation)

#  Display the conflicting regions where both helices and beta-strands are predicted.
def Display_The_Conflicting_Regions(helix_data, beta_strands_data, line_length):
    for position in range(0, len_of_sequence, line_length):
        start = position
        end = position + line_length
        protein = protein_sequence[start:end]
        print()
        print(protein)   # It prints the conflicting regions along with their corresponding secondary structure elements.
        helix_notation = Generate_Secondary_Structure(1, start, len(protein), helix_data)
        strand_notation = Generate_Secondary_Structure(2, start, len(protein), beta_strands_data)
        print(helix_notation)
        print(strand_notation)
        print()

# This function resolve conflicts between helices and beta-strands in the conflicting regions.
def resolveForHelixAndSheet(conflicting_sequence):
    helix_score=0
    strand_score=0
    answer=[]
    for i in conflicting_sequence:
        three_letter_code=One_letterCODE_to_Three_letterCODE[i]
        helix_score=helix_score+alpha[three_letter_code]
        strand_score=strand_score+beta[three_letter_code]
    # It calculates the total scores for helices and beta-strands in the conflicting sequence and assigns the dominant secondary structure element.
    if(helix_score>=strand_score):
        for i in range(len(conflicting_sequence)):
            answer.append(1)
    else:
        for i in range(len(conflicting_sequence)):
            answer.append(2)
    return answer


# This display the final assignment of secondary structure elements after resolving conflicts.
def Display_Final_Sequence_After_Solving_Conflicts(final_assignment,line_length):
    position=0
    while(position<len_of_sequence):
        print()
        start=position
        end=position+line_length
        protein=protein_sequence[start:end]
        print(protein)
        final=""
        for i in range(start,start+len(protein)):
            if(final_assignment[i]==0):
                final=final+" "
            elif(final_assignment[i]==1):
                final=final+"H"
            else:
                final=final+"S"
        print(final)
        position=end

# This resolve conflicts between helices and beta-strands in the conflicting regions.
def Resolving_Conflicts(helix_zero_one,strand_zero_one,line_length):
    final_assignment=[]
    # 0 for none , 1 for helix , 2 for strand
    index=0
    while(index<len(helix_zero_one)):
        helix_answer=helix_zero_one[index]
        strand_answer=strand_zero_one[index]
        if(helix_answer==0 and strand_answer==0):
            final_assignment.append(0)
            index=index+1
        elif ((helix_answer==0 and strand_answer==1) or (helix_answer==1 and strand_answer==0)):
            if(helix_answer==0):
                final_assignment.append(2)
            else:
                final_assignment.append(1)
            index=index+1
        else:
            conflicting_sequence=""       # iterates through the conflicting regions
            while(helix_answer==1 and strand_answer==1):
                conflicting_sequence=conflicting_sequence+protein_sequence[index]
                index=index+1
                if(index<len_of_sequence):
                    helix_answer=helix_zero_one[index]
                    strand_answer=strand_zero_one[index]
                else:
                    break
            resolved_answer=resolveForHelixAndSheet(conflicting_sequence)  # resolves the conflicts
            for i in resolved_answer:
                final_assignment.append(i)
    Display_Final_Sequence_After_Solving_Conflicts(final_assignment,line_length)  # generates the final assignment
    
# Finding nucleation sites for helices and beta-strands
helix_pointers=Finding_Nucleation_Sites(6,1,4)
strand_pointers=Finding_Nucleation_Sites(5,2,3)

# Extending nucleation sites for helices and beta-strands
final_helix=extendWindow(helix_pointers,4,1)
final_betaStrand=extendWindow(strand_pointers,4,2)

# Generating 0/1 lists for helices and beta-strands
helix_zero_one_list=Generate_0_1_List(final_helix)
strand_zero_one_list=Generate_0_1_List(final_betaStrand)

print()
print("---------------------------------------------------------------------------------------------------------------------------------------------------------")
print(" A) THE SEQUENCE REGIONS THAT ARE HELICAL IN NATURE .")
Display_The_Sequence_Regions(helix_zero_one_list,1,150)
print()
print("---------------------------------------------------------------------------------------------------------------------------------------------------------")

print(" B) THE SEQUENCE REGIONS THAT HAVE A TENDENCY TO FORM BETA STRANDS .")
Display_The_Sequence_Regions(strand_zero_one_list,2,150)
print()
print("---------------------------------------------------------------------------------------------------------------------------------------------------------")

print(" C) DISPLAYING THE CONFLICTING REGIONS .")
Display_The_Conflicting_Regions(helix_zero_one_list,strand_zero_one_list,150)
print("---------------------------------------------------------------------------------------------------------------------------------------------------------")

print(" # FINAL ASSIGNMENT OF SECONDARY STRUCTURAL ELEMENTS AFTER SOLVING CONFLICTS .")
Resolving_Conflicts(helix_zero_one_list,strand_zero_one_list,150)
print()
print("---------------------------------------------------------------------------------------------------------------------------------------------------------")
