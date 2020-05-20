import numpy as np

#test sequences
S1 = "ATATACTGG"
S2 = "CTAGGTGATATA"

# Create dynamic programming matrix
DPMatrix = np.ndarray(shape=(len(S1)+1,len(S2)+1), dtype=int)
DPMatrix.fill(0)

# Create matrix for back pointers
BackPtr = np.ndarray(shape=(len(S1)+1,len(S2)+1), dtype=int)
BackPtr.fill(-9)

#scoring matrix
# A:A, T:T, G:G, C:C get +1 any mismatch gets -1 and gap penalty gets -1
ScoreMat = {
'A': {'A': 1,'C':-1,'G':-1,'T':-1}, 
'C': {'A':-1,'C': 1,'G':-1,'T':-1}, 
'G': {'A':-1,'C':-1,'G': 1,'T':-1},
'T': {'A':-1,'C':-1,'G':-1,'T': 1},
}

indel_penalty = -1

# Filling up the DP matrix
for i in range(1,len(S1)+1):
    for j in range(1,len(S2)+1):
        char1 = S1[i-1] # current character at seq 1
        char2 = S2[j-1] # current character at seq 2

        # 0 for up, 1 for diag, 2 for left
        scores = [-999,-999,-999]
        # score diagonal
        scores[1] = DPMatrix[i-1,j-1] + ScoreMat[char1][char2]
        # score up: gap in sequence 2
        scores[2] = DPMatrix[i-1,j] + indel_penalty
       # score left: gap in sequence 1
        scores[0] = DPMatrix[i,j-1] + indel_penalty
        
        #if there is degeneracy, choose diagonal move
        if((scores[0]==scores[1] and scores[0]>scores[2] ) or (scores[2]==scores[1] and scores[2]>scores[0])):
            BestScore=scores[1]
        else:
            BestScore = max(scores)  # find the maximum score 
            
        DPMatrix[i,j]=BestScore
        for k in range(3):
            if scores[k] == BestScore:
                BackPtr[i,j] = k     #fill the position of maximum score in back pointer to mark up, diagonal or left

print("Dynamic programming matrix:")
print(DPMatrix)

print("\nBack pointers:")
print(BackPtr)

#Find maximum value in the last row and last column of DP matrix
i=len(S1)
maxrow=DPMatrix[i,0]
rowno=0
for j in range(len(S2)):
    if(maxrow<DPMatrix[i,j]):
        maxrow=DPMatrix[i,j]
        rowno=j

j=len(S2)
maxcol=DPMatrix[0,j]
colno=0
for i in range(len(S1)):
    if(maxcol<DPMatrix[i,j]):
        maxcol=DPMatrix[i,j]
        colno=i    
    
# Reading out the backtrace to get the optimal overlap alignment
align1 = ""
align2 = ""

#If the maximum value is in the last column of DP Matrix
if(maxcol>=maxrow): #prefix of S1 aligns with suffix of S2
    
    i=colno
    j=len(S2)
    str1=""
    str2=""
    
    while i>0 or j>0:
        if BackPtr[i,j] == 0: # gap in sequence 1
            align1 += "-"
            align2 += S2[j-1]
            j -= 1
            if(i==0 or j==0): break
        if BackPtr[i,j] == 1: # diagonal
            align1 += S1[i-1]
            align2 += S2[j-1]
            i -= 1
            j -= 1
            if(i==0 or j==0): break
        if BackPtr[i,j] == 2: # gap in equence 2
            align1 += S1[i-1]
            align2 += "-"
            i -= 1
            if(i==0 or j==0): break
   
    #Display the optimal score of alignment    
    print("\n")
    print(align1)
    print(align2)
    print("Optimal score of alignment= " + str(maxcol))
    
    #Adding remaining part of sequences and/or periods(.) to the aligned sequences
    k=len(S1)-1
    while k>0:
        if(align1[len(align1)-1]!=S1[k]):
            str1+=S1[k]
            k-=1
        else:
            break
    align1=align1[::-1]
    str11=str1[::-1]
    align1+=str11
    if(S2.find(align2)):
        pos2=S2.find(align2)
        str2+=S2[0:pos2-1]
    str22=str2[::-1]
    align2+=str22
    align2=align2[::-1]
    align2+='.'*len(str11)
    align1='.'*len(str22)+align1
    
    
#If the maximum value is in the last row of DP Matrix
if(maxcol<maxrow):  #suffix of S1 aligns with prefix of S2 
    i=len(S1)
    j=rowno   
    str1=""
    str2=""
   
    while i>1 or j>1:
        if BackPtr[i,j] == 0: # gap in sequence 1
            align1 += "-"
            align2 += S2[j-1]
            j -= 1
            if(i==0 or j==0): break
        if BackPtr[i,j] == 1: # diagonal
            align1 += S1[i-1]
            align2 += S2[j-1]
            i -= 1       
            j -= 1
            if(i==0 or j==0): break
        if BackPtr[i,j] == 2: # gap in equence 2
            align1 += S1[i-1]
            align2 += "-"
            i -= 1 
            if(i==0 or j==0): break
    
    #Display the optimal score of alignment        
    print("\n")       
    print(align1[::-1])
    print(align2[::-1])
    print("Optimal score of alignment= " + str(maxrow))        
    
    #Adding remaining part of sequences and/or periods(.) to the aligned sequences    
    for i in range(len(S1)):
        if(align1[len(align1)-1]!=S1[i]):
            str1+=S1[i]
        else:
            break
    align1+=str1
    align2+='.'*len(str1)
    align1=align1[::-1]
    align2=align2[::-1] 
    S22=S2[::-1]
    for k in range(len(S1)):
        if(align2[len(align2)-1]!=S2[k]):
            str2+=S2[k]
        else:
            break
    str22=str2[::-1]
    align2+=str22

print("\nOptimal overlap alignment between the sequences:")
print(align1)
print(align2)