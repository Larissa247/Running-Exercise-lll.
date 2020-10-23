# -*- coding: utf-8 -*-
"""
The program should be run as MSA.py fasta_file weight_parameters output_file

What is meant by weight parameters? 
--> gives user the freedom to determine which score is assigned to a match/no-match?
        is that what is meant? YES
        should weight parameters be a txt. file? YES
        
Procedure:
    1. Input file GeneticData.txt is manually edited by adding a '>' sign in front of each name
    2. Iterate through GeneticData.txt file and put mitochondrial DNA and Y chromosomal 
        DNA sequences in two different output files (fasta format) including their respective names
        
    3. Create two dictionaries:
        - mtDNA_dict which contains names as key and mtDNA as value
        - YchrDNA_dict which contains names as key and Ychromosomal DNA as value
    4. Read weight_parameter input file in which the user could define the penalty scores
        - Create new list for penalty scores (called parameter_list)
        - Iterate through weight_parameter file
        - and append each score (as integer) to parameter_list
        - assign variables (that can be called in function) to each item of the parameter_list
    5. Create a function (score_one_pair) that defines the scoring system for pairwise alignments:
        - score for gap and gap/or "?" and "?" = 0
        - score for gap to nucleotide = -1
        - score for same nucleotide (match) = 1
        - score for A-G or C-T transition = -1
        - score for transversion = -2
        - score for occurance of "?" = -1
        
    6. Define a function (score_two_ids) that scores to samples/Ids with each other:
        - the function makes use of the dictionary - obtains DNA sequence by calling value of dict.
        - zip() is used two pair each name/id with all the others, but only once! 
        
            (zip() returns a zip object, which is an iterator of tuples where the first item 
             in each passed iterator is paired together, and then the second item in each 
             passed iterator paired together ect.)
            
        - introduce a score and set it to zero (score =  0)
        - call the function from above (score_one_pair) within the new function score_two_ids,
            and calculate the score, finally return it
            
    7. Create a function that calculates the identity score of two sequences:
        - Note: big gaps should be ignored
        - 
    
    8. Create a new output file (MAS_results) that contains 


Created on Fri Oct 23 10:18:53 2020

@author: Larissa




"""

#Write a python code that will do multiple sequence alignment (MSA) of DNA sequences


#open GeneticData.txt file ('r' mode)
#extract mtDNA (with respective name) and export to new output file (fasta format) e.g. mtDNA.fna
#repeat procedure with Y chromosomal DNA e.g. export to YchrDNA.fna
    

with open('GeneticData.txt','r') as file,\
    open('mtDNA.fna','w') as mtDNA,\
        open('YchrDNA.fna','w') as YchrDNA:
            
            for line in file:                                           #iterate through lines in file
                if line.startswith('>'):                                #if line starts with '>' sign:
                    name = line.strip()                                 #take that line, remove leading and trailing 
                                                                        #characters with strip(), assign it variable "name"
                
                if 'mtDNA' in line:                                     #if mtDNA occurs in line:
                    sequence1 = next(file)                              #jump to next line in file (which is mtDNA seq)
                    sequence1 = sequence1.strip()                       #again remove leading/trailing characters
                    
                    mtDNA.write('{}\n{}\n'.format(name,sequence1))      #append name,newline,mtDNA-seq to ouput file 1
                
                if 'Y chromosome' in line:                              #if Y chromosome occurs in line:
                    sequence2 = next(file)                              #jump to next line in file
                    sequence2 = sequence2.strip()                       #apply strip
                    
                    YchrDNA.write('{}\n{}\n'.format(name,sequence2))    #append name,newline,YchrDNA-seq to output file 2




#Create dictionaries and fill them with input from mtDNA and YchrDNA fasta files

#mtDNA_dict.clear()
#YchrDNA_dict.clear()
    
mtDNA_dict = dict()
YchrDNA_dict = dict()

with open('mtDNA.fna','r') as mtDNA,\
    open('YchrDNA.fna','r') as YchrDNA:
        
        for line in mtDNA:
            if line.startswith('>'):                    #if line startswith '>':
                NameIsKey = line[1:].strip()            #remove '>', strip \n away and define that line as key
            else:                                       #else:
                SeqIsValue = line.strip()               #strip \n away and define that line as value
             
                mtDNA_dict[NameIsKey] = SeqIsValue      #add name as key and sequence as value to dictionary


        for line in YchrDNA:
            if line.startswith('>'):                    #if line startswith '>':
                NameIsKey = line[1:].strip()            #remove '>', strip \n away and define that line as key
            else:                                       #else:
                SeqIsValue = line.strip()               #strip \n away and define that line as value
             
                YchrDNA_dict[NameIsKey] = SeqIsValue    #add name as key and sequence as value to dictionary          



#works! But probably "nicer" to create a function for this instead of repeating code


#Reading the users input for weight parameters (penalties) and assigning variables to it

parameter_list = []                                     #Create parameter list
#parameter_list.clear()

with open('weight_parameters.txt','r') as wp:           #open weight_parameter file in read mode
    for line in wp:                                     #iterate through lines in file
        if line.startswith('Score'):                    #if line starts with Score:
            users_input = next(wp)                      #jump to next line and define it as users_input
            users_input = int(users_input.strip())      #strip and convert string to integer
            parameter_list.append(users_input)          #append users_input to parameter_list
                        
wp.close()                                              #close file

#Assign variables to list elements
GapGap = parameter_list[0]
GapOther = parameter_list[1]
UnknownUnknown = parameter_list[2]
UnknownNt = parameter_list[3]
match = parameter_list[4]
transition = parameter_list[5]
transversion = parameter_list[6]


#Introducing a scoring systems 

#Introduce scoring system for one-to-one position:
#return values are defined by user (see above)


def score_one_pair(x, y):               #defines function for scoring two characters, e.g. in a tuple 
    if x == '-' and y == '-':           #both x and y are -
        return GapGap                   #score = 0
    
    if x == '-' or y == '-':            #x or y are -. Both could be, but loop will break and return score of previous line
        return GapOther                 #score = -1
    
    if x == '?' and y == '?':           #same procedure for '?' as '-'
        return UnknownUnknown           #score = 0
    
    if x == '?' or y == '?':
        return UnknownNt                #score = -1
    
    if x == y:                          #x,y are the same character, but neither can be -
        return match                    #score = 1
    
    
    # A-G Transition
    if (x == "A" and y == "G") or (x == "G" and y == "A"):
        return transition                                           #score = -1
    
    # C-T Transition
    if (x == "C" and y == "T") or (x == "T" and y == "C"):
        return transition                                           #score = -1
    
    # transversions A <-> CT, G <-> CT 
    transversions = {"A": ["C", "T"], "G": ["C", "T"], "C": ["A", "G"], "T": ["A", "G"]}
    #tranversions is now a dictionary where one key corresponds to two values
    if y in transversions[x]: #if y is the value of x in transversions 
        return transversion                 #score = -2

    print(y)
    if x in transversions[y]: #if x is the value of the transversions 
        return transversion                 #score = -2



#Define function for scoring 2 sequences (stored in dictionary) -
#puting first seq in x and second seq in y:

def score_two_names(random_dict,firstName,secondName):          
    firstDna = random_dict[firstName]                       #first DNA-seq is defined as value of respective key in dictionary 
    secondDna = random_dict[secondName]                     #second DNA-seq is defined as value of respective key in dictionary
    both = zip(firstDna, secondDna)                         #both will be tuples of the integers; if the sequences would be of different lenght, 
                                                            #it would align them in a tuple until the shortest sequence is exhausted

    score = 0                                               #introduce score and set to 0
    for [x,y] in both:                                      #for one character in the tuples (the sequences) at the time, until the end of the tuple
        score = score + score_one_pair(x,y)                 #uses scoring of one-pair-function, then adds that to the existing score
    return score


#Calculate identity score
#Define function for scoring two characters:

def identical_nts(x, y):                            #x represents character 1, y represents character 2
    if x == '-' and y == '-':                       #if x and y are both "-"
       return 0                                     #return 0
    if x == '?' and y == '?':                       #if x and y are both "?"
        return 0                                    #return 0
    if x == y:                                      #if x is the same as y (identical/match!)
       return 1                                     #return 1
    if not x == y:                                  #if x and y are not the same
       return 0                                     #return 0


#Define function that calculates identity score for two DNA-sequences:

def score_two_samples(random_dict,firstSample,secondSample):
    seq1 = random_dict[firstSample]                                 #seq1 = value of key (firstSample) in dictionary (random_dict) 
    seq2 = random_dict[secondSample]                                #seq2 = value of key (secondSample) in dictionary (random_dict)
    both = zip(seq1,seq2)                                           #perform zip() function on both sequences
    
    a = ""                                                          #define a as empty string
    b = ""                                                          #define b as empty string
    identical = 0                                                   #introduce counter "identical", set it to 0
    for [x, y] in both:                                             #for each character(x,y) in the sequences
        identical = identical + identical_nts(x, y)                     #call function that determines amount of identical nucleotides and add to counter
        if not (x == "-" and y == "-") or (x == "?" and y == "?"):      #if both positions are not gaps (-) or unknown (?): 
            a += x                                                      #append character(x) to empty string a
            b += y                                                      #append character(y) to empty string b
    length = len(a)                                                 #define length as length of string a
    identityScore = 100* (identical/length)                         #calculate percentage
    identityScore = identityScore
    return identityScore



#Write output file and iterate through key is dictionary to obtain scores:

f = open('MSA_results.txt','w')
f.write("")
f.close()

with open('MSA_results.txt', "a") as file:
    
    def give_all_scores(choose_dict):
        new_list = list(choose_dict)
        identity_score_list = list()
        score_list = list()
        n = len(new_list)
        
        for i in range(n):
            for j in range(i+1,n):
                firstName = new_list[i]
                secondName = new_list[j]
                firstSeq = choose_dict[firstName]
                secondSeq = choose_dict[secondName]
                identity_score = score_two_samples(choose_dict,firstName,secondName)
                score = score_two_names(choose_dict,firstName,secondName)
                
                identity_score = str("{:.2f}".format(identity_score))
                score = str(score)
                
                file.write(firstName + "\t" + secondName + "\t" + identity_score + "%"+"\t" + score + "\n")
                
                identity_score_list.append(identity_score)
                score_list.append(score)
            
            return identity_score_list, score_list 
            

mtDNA_scores = give_all_scores(mtDNA_dict)
YchrDNA_scores = give_all_scores(YchrDNA_dict)



