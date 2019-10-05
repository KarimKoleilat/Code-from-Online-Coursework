__author__ = 'Karim'
"""
Sequence alignment problem
"""

DESKTOP = True

import math
import random
import urllib2
import math
if DESKTOP:
    import matplotlib.pyplot as plt
    #import alg_project4_solution as student
else:
    import simpleplot
    import userXX_XXXXXXX as student


# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"



###############################################
# provided code

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urllib2.urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict








def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urllib2.urlopen(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq





def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urllib2.urlopen(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print "Loaded a dictionary with", len(word_list), "words"
    return word_list



def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    function takes an alphabet and scoring values and returns a scoring matrix
    in the form of a dictionary of dictionaries.
    """
    copy_alph = alphabet.copy()
    copy_alph.add('-')
    matrix_score = {}
    for char in copy_alph:
        matrix_score[char]={}
        for char2 in copy_alph:
            if char == '-' or char2 == '-':
                matrix_score[char][char2] = dash_score
            elif char == char2:
                matrix_score[char][char2] = diag_score
            else:
                matrix_score[char][char2] = off_diag_score

    return matrix_score



def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    returns the alignment matrix of two sequences with option to do local or
    global pairwise alignment.
    """
    len_x = len(seq_x)+1
    len_y = len(seq_y)+1
    align_matrix = [[0 for dummy_idy in range(len_y)] for dummy_idx in range(len_x)]


    for idx in range(1,len_x):

        if global_flag :
            align_matrix[idx][0] = scoring_matrix[seq_x[idx-1]]['-'] + align_matrix[idx-1][0]

    for idy in range(1,len_y):

        if  global_flag :
            align_matrix[0][idy] = scoring_matrix['-'][seq_y[idy-1]] + align_matrix[0][idy-1]


    for idx in range(1,len_x):
        for idy in range(1,len_y):
            val_1 = align_matrix[idx-1][idy-1]+scoring_matrix[seq_x[idx-1]][seq_y[idy-1]]
            val_2 = align_matrix[idx-1][idy] + scoring_matrix[seq_x[idx-1]]['-']
            val_3 = align_matrix[idx][idy-1] + scoring_matrix['-'][seq_y[idy-1]]
            val_to_add = max(val_1,val_2,val_3)

            if (not global_flag) and val_to_add<0:
                #val_1 = align_matrix[idx-1][idy-1]
                #val_2 = align_matrix[idx-1][idy]
                #val_3 = align_matrix[idx][idy-1]
                val_to_add = 0

            align_matrix[idx][idy] = val_to_add

    if not global_flag:
        for idx in range(len_x):
            for idy in range(len_y):
                if align_matrix[idx][idy] <0:
                    align_matrix[idx][idy] =0

    """
    for idx in range(len_x):
        print align_matrix[idx]
    """
    return align_matrix

def compute_global_alignment(seq_x, seq_y, scoring_matrix, align_matrix):
    """
    Computes the global alignment of two sequences
    seq_x and seq_y
    """
    ind_x = len(seq_x)
    ind_y = len(seq_y)
    str_x = ""
    str_y = ""
    while ind_x >0 and ind_y >0:

        if align_matrix[ind_x][ind_y] == (align_matrix[ind_x-1][ind_y-1] + scoring_matrix[seq_x[ind_x-1]][seq_y[ind_y-1]]):
            str_x = seq_x[ind_x-1] + str_x
            str_y = seq_y[ind_y-1] + str_y
            ind_x = ind_x -1
            ind_y = ind_y -1
        elif align_matrix[ind_x][ind_y] == (align_matrix[ind_x-1][ind_y] + scoring_matrix[seq_x[ind_x -1 ]]['-']):
            str_x = seq_x[ind_x-1] + str_x
            str_y = '-' +str_y
            ind_x = ind_x -1
        else:

            str_y = seq_y[ind_y-1] +str_y
            str_x = '-' + str_x
            ind_y = ind_y -1

    while ind_x >0:

        str_x = seq_x[ind_x-1] + str_x

        str_y = '-' +str_y
        ind_x = ind_x -1

    while ind_y>0:
        str_x = '-' + str_x
        str_y = seq_y[ind_y-1] + str_y
        ind_y = ind_y - 1

    score = 0

    for idx in range(len(str_x)):
        score += scoring_matrix[str_x[idx]][str_y[idx]]

    return (score , str_x,str_y )

def max_indexes(align_matrix):
    """
    takes a 2d array and returns the indexes of the max value
    in the 2d array along with the value in the form of a tuple
    (value, x , y)

    """
    rows = len(align_matrix)
    cols = len(align_matrix[0])
    max_val = 0
    max_idx = 0
    max_idy = 0
    for idx in range(rows):
        for idy in range(cols):
            if align_matrix[idx][idy]>max_val:
                max_val = align_matrix[idx][idy]
                max_idx = idx
                max_idy = idy
    return (max_val,max_idx,max_idy)




def compute_local_alignment(seq_x, seq_y, scoring_matrix, align_matrix):
    """
    returns local alignment of two sequences
    """

    max_ref = max_indexes(align_matrix)
    ind_x = max_ref[1]
    ind_y = max_ref[2]
    str_x = ""
    str_y = ""
    while align_matrix[ind_x][ind_y]>0:

        if align_matrix[ind_x][ind_y] == (align_matrix[ind_x-1][ind_y-1] + scoring_matrix[seq_x[ind_x-1]][seq_y[ind_y-1]]):
            str_x = seq_x[ind_x-1] + str_x
            str_y = seq_y[ind_y-1] + str_y
            ind_x = ind_x -1
            ind_y = ind_y -1
        elif align_matrix[ind_x][ind_y] == (align_matrix[ind_x-1][ind_y] + scoring_matrix[seq_x[ind_x -1 ]]['-']):
            str_x = seq_x[ind_x-1] + str_x
            str_y = '-' +str_y
            ind_x = ind_x -1
        else:

            str_y = seq_y[ind_y-1] +str_y
            str_x = '-' + str_x
            ind_y = ind_y -1


    score = 0

    for idx in range(len(str_x)):
        score += scoring_matrix[str_x[idx]][str_y[idx]]

    return (score , str_x,str_y )

def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    dict_to_return = {}

    for idx in range(num_trials):
        y_copy = [char for char in seq_y]
        random.shuffle(y_copy)
        y_copy = ''.join(y_copy)
        trial_align_matrix = compute_alignment_matrix(seq_x,y_copy,scoring_matrix,False)
        trial_score = compute_local_alignment(seq_x, y_copy, scoring_matrix, trial_align_matrix)
        if dict_to_return.has_key(trial_score[0]):
            dict_to_return[trial_score[0]] = dict_to_return[trial_score[0]] + 1
        else:
            dict_to_return[trial_score[0]] = 1

    return dict_to_return

def compute_edit_distance(seq_x,seq_y,alphabet):
    len_x = len(seq_x)
    len_y = len(seq_y)

    diag_score = 2
    off_diag_score = 1
    dash_score = 0

    score_m = build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score)
    alignment_matrix = compute_alignment_matrix(seq_x, seq_y, score_m, True)
    global_a = compute_global_alignment(seq_x,seq_y,score_m,alignment_matrix)

    edit_distance = len_x + len_y - global_a[0]

    return edit_distance


def check_spelling(checked_word, dist, word_list):
    list_len = len(word_list)
    alphabet = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'}
    list_to_ret = []
    for word in word_list:
        if compute_edit_distance(checked_word,word,alphabet)<= dist:
            list_to_ret.append(word)

    return list_to_ret


str_a = "ACGTCATTTA"
str_b = "ACGTCATTAA"

#print "this is the edit_distance: " + str(compute_edit_distance(str_a,str_b))
file = open('assets_scrabble_words3.txt')
word_list = file.read()
word_list = word_list.splitlines()
#print word_list

print check_spelling("firefly", 2, word_list)




file.close()




scoring_PAM50 = read_scoring_matrix(PAM50_URL)
consensus_pax = read_protein(CONSENSUS_PAX_URL)
human_eyeless =  read_protein(HUMAN_EYELESS_URL)
fruitfly_eyeless = read_protein(FRUITFLY_EYELESS_URL)

aa_alphabet = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','-'}
#prot_alignment_matrix = compute_alignment_matrix(human_eyeless, fruitfly_eyeless, scoring_PAM50 , False)
#local_align_hum_fly = compute_local_alignment(human_eyeless, fruitfly_eyeless, scoring_PAM50, prot_alignment_matrix)

#local_hum = local_align_hum_fly[1].replace('-','')
#local_fruitfly = local_align_hum_fly[2].replace('-','')

#print "score is: " + str(local_align_hum_fly[0])
#hum_concensus_align = compute_alignment_matrix(local_hum, consensus_pax, scoring_PAM50, True)
#fruitfly_concensus_align = compute_alignment_matrix(local_fruitfly, consensus_pax, scoring_PAM50, True)

#global_hum_concensus = compute_global_alignment(local_hum, consensus_pax, scoring_PAM50, hum_concensus_align )
#len_hum = len(global_hum_concensus[1])
#count_hum = 0
#for idx in range(len_hum):
#    if global_hum_concensus[1][idx] == global_hum_concensus[2][idx]:
#        count_hum +=1

#hum_pax_percentage = (float(count_hum)/len_hum)*100
#print "human percentage is: " + str(hum_pax_percentage)


#global_fruitfly_concensus = compute_global_alignment(local_fruitfly, consensus_pax, scoring_PAM50, fruitfly_concensus_align )
#len_fruitf = len(global_fruitfly_concensus[1])
#count_fruitf = 0
#for idx in range(len_fruitf):
#    if global_fruitfly_concensus[1][idx] == global_fruitfly_concensus[2][idx]:
#        count_fruitf +=1

#fruitf_pax_percentage = (float(count_fruitf)/len_fruitf)*100
#print "fruitfly percentage is: " + str(fruitf_pax_percentage)

#print "probability is: " + str(float(23**(-132)))
#print "human local alignment is: " + str(local_hum)
#print "fruitfly local alignment is: " + str(local_fruitfly)

"""
null_distribution = generate_null_distribution(human_eyeless,fruitfly_eyeless,scoring_PAM50,1000)
print null_distribution.keys()
print null_distribution
x_vals = null_distribution.keys()
y_vals = []
for xval in x_vals:
    val_to_add = float(null_distribution[xval])/1000
    y_vals.append(val_to_add )
plt.ylabel( "fraction of the frequency of a given alignment score" )
plt.xlabel( "alignment score of hum_local_alignment with shuffled fruitfly_local_alignment" )
plt.title("Normalized distribution of the frequencies of each alignment score of hum_local and shuffled fruitfly_local")
plt.bar(x_vals,y_vals,0.83,color="blue")
plt.show()

distribution_dict = {36: 1, 38: 1, 39: 4, 40: 6, 41: 5, 42: 26, 43: 31, 44: 53, 45: 50, 46: 57, 47: 66, 48: 48, 49: 64, 50: 74, 51: 73, 52: 54, 53: 60, 54: 39, 55: 52, 56: 30, 57: 27, 58: 33, 59: 17, 60: 27, 61: 26, 62: 10, 63: 7, 64: 9, 65: 11, 66: 7, 67: 10, 68: 4, 69: 5, 70: 2, 73: 1, 74: 4, 75: 2, 78: 1, 80: 1, 81: 1, 88: 1}
scores = distribution_dict.keys()
mean_scores = 0

for score in scores:
    mean_scores+=score


mean_scores = float(mean_scores)/1000
dev_sum = 0
for score in scores:
    dev_sum+=(score-mean_scores)**2

std_dev = (0.001*dev_sum)**.5
local_alignment_score = 875
z_score = float(local_alignment_score-mean_scores)/std_dev

print "The mean is: " +str(mean_scores)
print "The standard deviation is: " +str(std_dev)
print "The z-score for the alignment score is: " +str(z_score)









seq_1 = 'TGCTTC'
seq_2 =  'GGACTGCTTCTGG'

seq_3 =  "GAACTGCTTCTGG"

score_m1 = build_scoring_matrix(set(['A','C','G','T']), 5, 2, -4)
score_m2 = {'A': {'A': 2, 'C': 1, '-': 0, 'T': 1, 'G': 1},
            'C': {'A': 1, 'C': 2, '-': 0, 'T': 1, 'G': 1},
            '-': {'A': 0, 'C': 0, '-': 0, 'T': 0, 'G': 0},
            'T': {'A': 1, 'C': 1, '-': 0, 'T': 2, 'G': 1},
            'G': {'A': 1, 'C': 1, '-': 0, 'T': 1, 'G': 2}}

align_matrix1 = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                [0, 1, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                [0, 1, 2, 3, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                [0, 1, 2, 4, 4, 6, 7, 7, 7, 7, 7, 7, 7, 7],
                [0, 1, 2, 4, 6, 6, 7, 9, 9, 9, 9, 9, 9, 9],
                [0, 1, 2, 4, 6, 8, 8, 9, 11, 11, 11, 11, 11, 11]]

result = compute_local_alignment(seq_1, seq_2, score_m2, align_matrix1)

print result
#print seq_2




reference = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 2, 6, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0], [0, 7, 5, 3, 0, 0, 0, 0, 0, 7, 2, 0, 0, 0, 0, 0], [0, 3, 14, 9, 4, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0], [0, 0, 9, 20, 15, 10, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 4, 15, 29, 24, 19, 14, 9, 4, 0, 0, 0, 0, 0, 0], [0, 0, 0, 10, 24, 37, 32, 27, 22, 17, 12, 7, 2, 0, 0, 3], [0, 0, 0, 5, 19, 32, 46, 41, 36, 31, 26, 21, 16, 11, 6, 1], [0, 0, 0, 0, 14, 27, 41, 52, 47, 42, 37, 32, 27, 22, 17, 12], [0, 0, 0, 0, 9, 22, 36, 47, 46, 44, 39, 34, 29, 24, 19, 17], [0, 0, 0, 0, 4, 17, 31, 42, 41, 42, 38, 33, 28, 24, 19, 18], [0, 0, 2, 0, 0, 12, 26, 37, 49, 44, 39, 39, 34, 29, 24, 19], [0, 2, 0, 0, 0, 7, 21, 32, 44, 57, 52, 47, 42, 37, 32, 27], [0, 0, 0, 0, 0, 2, 16, 27, 39, 52, 65, 60, 55, 50, 45, 40], [0, 0, 0, 0, 0, 0, 11, 22, 34, 47, 60, 71, 66, 61, 56, 51], [0, 0, 0, 0, 0, 0, 6, 17, 29, 42, 55, 66, 79, 74, 69, 64], [0, 0, 0, 0, 0, 0, 1, 12, 24, 37, 50, 61, 74, 85, 80, 75], [0, 0, 0, 0, 0, 0, 0, 7, 19, 32, 45, 56, 69, 80, 98, 93], [0, 0, 0, 0, 0, 3, 0, 2, 14, 27, 40, 51, 64, 75, 93, 105], [0, 0, 0, 0, 3, 0, 1, 0, 9, 22, 35, 46, 59, 70, 88, 100], [0, 0, 0, 0, 0, 0, 0, 0, 4, 17, 30, 41, 54, 65, 83, 95], [0, 6, 2, 0, 0, 0, 0, 0, 0, 12, 25, 36, 49, 60, 78, 90]]
str1 = ""
for idy in range(len(reference)):
    if reference[idy] != result[idy]:
        str1 += str(idy) +"  "
print
print str1
print

for idx in reference:
    print idx


len_1 = len(seq_1)
len_2 = len(seq_2)

for idx in range(len_1):
    for idy in range(len_2):
        if reference[idx][idy] != result[idx][idy]:
            print idx, idy
            print "found the problem"
            print

if reference == result:
    print "it's working?!"
else:
    print "something's wrong"
"""


