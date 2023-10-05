import sys
import copy, enum

E = 0.02 # 엡실론

class Nucleotide(enum.Enum):
    A = 0
    C = 1
    G = 2
    T = 3

def open_and_parse(input_file: str):
    f = open(input_file, 'r')

    offset_hashmap = dict()
    titles = list()
    max_len = 0
    seq = ''
    while True:
        line = f.readline().rstrip()

        if not line:
            if seq:
                length = len(seq)
                max_len = max(length, max_len)

            f.close()
            break
        
        elif line[0] == '>' and seq:
            length = len(seq)
            max_len = max(length, max_len)

            seq = ''

            title = line[1:]
            offset_hashmap[title] = f.tell()
            titles.append(title)
        elif line[0] != '>':
            seq += line
        elif line[0] == '>':
            title = line[1:]
            offset_hashmap[title] = f.tell()
            titles.append(title)
        

    return titles, offset_hashmap, max_len



def process_reverse(forward_strand: str):
    reversed_strand = ""
    for c in forward_strand:
        if c == 'A':
            reversed_strand += 'T'
        elif c == 'T':
            reversed_strand += 'A'
        elif c == 'C':
            reversed_strand += 'G'
        elif c == 'G':
            reversed_strand += 'C'
    
    reversed_strand = reversed_strand[::-1]
    return reversed_strand

def update_hidden(hidden: list, profile: list, titles: list, offset_hashmap: dict, input_file: str, motif_len: int):
    f = open(input_file, 'r')

    for i in range(len(titles)):
        seq_title = titles[i]
        f.seek(offset_hashmap[seq_title])

        nucleotide_seq = ''
        while True:
            line = f.readline()
            if not line or line[0] == '>':
                break
            nucleotide_seq += line.rstrip()
        
        for j in range(len(nucleotide_seq)-motif_len+1):
            probability = 1
            for k in range(motif_len):
                if nucleotide_seq[k+j] == 'A':
                    probability *= profile[Nucleotide.A.value][k]
                elif nucleotide_seq[k+j] == 'C':
                    probability *= profile[Nucleotide.C.value][k]
                elif nucleotide_seq[k+j] == 'G':
                    probability *= profile[Nucleotide.G.value][k]
                elif nucleotide_seq[k+j] == 'T':
                    probability *= profile[Nucleotide.T.value][k]
                else:
                    print('[Error]')
                    exit(-1)
            hidden[i][j] = probability

        col_num = len(nucleotide_seq)-motif_len+1
        reversed_strand = process_reverse(nucleotide_seq)
        for j in range(col_num):
            probability = 1
            for k in range(motif_len):
                if reversed_strand[k+j] == 'A':
                    probability *= profile[Nucleotide.A.value][k]
                elif reversed_strand[k+j] == 'C':
                    probability *= profile[Nucleotide.C.value][k]
                elif reversed_strand[k+j] == 'G':
                    probability *= profile[Nucleotide.G.value][k]
                elif reversed_strand[k+j] == 'T':
                    probability *= profile[Nucleotide.T.value][k]
                else:
                    print('[Error]')
                    exit(-1)
            hidden[i][j+col_num] = probability
        
        total = sum(hidden[i])
        for j in range(len(hidden[0])):
            hidden[i][j] /= total


    return hidden

def update_profile(hidden: list, profile: list, titles: list, offset_hashmap: dict, input_file: str, motif_len: int, eps:float):
    # 프로필에 엡실론 값 더해주기
    for row in profile:
        for j in range(len(row)):
            row[j] = eps

    for i in range(len(titles)):
        seq_title = titles[i]
        f = open(input_file, 'r')
        f.seek(offset_hashmap[seq_title])

        nucleotide_seq = str()
        while True:
            line = f.readline()
            if not line or line[0] == '>':
                f.close()
                break
            nucleotide_seq += line.rstrip()
        for j in range(len(nucleotide_seq)-motif_len+1):
            for k in range(motif_len):
                if nucleotide_seq[j+k] == 'A':
                    profile[Nucleotide.A.value][k] += hidden[i][j]
                elif nucleotide_seq[j+k] == 'C':
                    profile[Nucleotide.C.value][k] += hidden[i][j]
                elif nucleotide_seq[j+k] == 'G':
                    profile[Nucleotide.G.value][k] += hidden[i][j]
                elif nucleotide_seq[j+k] == 'T':
                    profile[Nucleotide.T.value][k] += hidden[i][j]
                else:
                    print('[Error]')
                    exit(-1)
        
        reversed_strand = nucleotide_seq[::-1]
        col_num = len(nucleotide_seq)-motif_len+1
        for j in range(len(nucleotide_seq)-motif_len+1):
            for k in range(motif_len):
                if reversed_strand[j+k] == 'A':
                    profile[Nucleotide.A.value][k] += hidden[i][j+col_num]
                elif reversed_strand[j+k] == 'C':
                    profile[Nucleotide.C.value][k] += hidden[i][j+col_num]
                elif reversed_strand[j+k] == 'G':
                    profile[Nucleotide.G.value][k] += hidden[i][j+col_num]
                elif reversed_strand[j+k] == 'T':
                    profile[Nucleotide.T.value][k] += hidden[i][j+col_num]
                else:
                    print('[Error]')
                    exit(-1)
    
    # column 다 더해서 나누어주기
    for col in range(motif_len):
        total = 0
        
        for row in range(4):
            total += profile[row][col]

        for i in range(4):
            profile[i][col] /= total
    f.close()

    return profile

# 현재 프로필과 과거 프로필이 같은지를 확인
def is_same_profile(cur_profile, prev_profile):
    for i in range(len(cur_profile)):
        for j in range(len(cur_profile[i])):
            if cur_profile[i][j] != prev_profile[i][j]:
                return False
    
    return True


def init_matrix(profile:list, hidden:list, titles: list, offset_hashmap: dict, input_file:str,  motif_len, eps: float):
    # 프로필에 엡실론 값 더해주기
    for row in profile:
        for j in range(len(row)):
            row[j] = eps
    
    for i in range(len(titles)):
        seq_title = titles[i]
        f = open(input_file, 'r')
        f.seek(offset_hashmap[seq_title])

        nucleotide_seq = str()
        while True:
            line = f.readline()
            if not line or line[0] == '>':
                f.close()
                break
            nucleotide_seq += line.rstrip()
        
        for j in range(len(nucleotide_seq)-motif_len+1):
            for k in range(motif_len):
                if nucleotide_seq[j+k] == 'A':
                    profile[0][k] += 1
                elif nucleotide_seq[j+k] == 'C':
                    profile[1][k] += 1
                elif nucleotide_seq[j+k] == 'G':
                    profile[2][k] += 1
                elif nucleotide_seq[j+k] == 'T':
                    profile[3][k] += 1
                else:
                    print('[Error]')
                    exit(-1)
        
        col_num = len(nucleotide_seq)-motif_len+1
        reversed_strand = process_reverse(nucleotide_seq)
        for j in range(col_num):
            for k in range(motif_len):
                if reversed_strand[j+k] == 'A':
                    profile[0][k] += 1
                elif reversed_strand[j+k] == 'C':
                    profile[1][k] += 1
                elif reversed_strand[j+k] == 'G':
                    profile[2][k] += 1
                elif reversed_strand[j+k] == 'T':
                    profile[3][k] += 1
                else:
                    print('[Error]')
                    exit(-1)


    # column 다 더해서 나누어주기
    for col in range(motif_len):
        total = 0
        
        for row in range(4):
            total += profile[row][col]

        for i in range(4):
            profile[i][col] /= total
    # for i in range(4):
    #     profile[i][0] = 1
    
    # hidden matrix 만들기
    hidden = update_hidden(hidden, profile, titles, offset_hashmap, input_file, motif_len)

    return profile, hidden


def find_motifs(hidden, input_file, titles, offset_hashmap, motif_length):
    motifs_array = []
    for i in range(len(titles)):
        max_idx, max_prob = 0, 0
        for j in range(len(hidden[i])):
            if max_prob < hidden[i][j]:
                max_idx = j
                max_prob = hidden[i][j]
        
        f = open(input_file, 'r')
        f.seek(offset_hashmap[titles[i]])
        seq = ''
        while True:
            line = f.readline().rstrip()
            if not line or line[0] == '>':
                f.close()
                break
            seq += line
        col_num = len(seq)
        if max_idx < col_num:
            motif = seq[max_idx: max_idx+motif_length]
            # if not motif:
            #     print("line 278: ", max_idx, motif_length, col_num)
            #     exit(1)
        else:
            rev = process_reverse(seq)
            motif = rev[max_idx-col_num: max_idx+motif_length-col_num]
            # if not motif:
            #     print("line 284: ", max_idx, motif_length, col_num)
            #     exit(1)

        motifs_array.append(motif)
    return motifs_array
            



if __name__=='__main__':
    if len(sys.argv) < 3:
        print("[ERROR]: Insufficient args")
        sys.exit()
    
    input_file = sys.argv[1]
    k = int(sys.argv[2]) # 보통 11, 12 mer
    titles, offset_hashmap, n = open_and_parse(input_file) # n은 max len
    
    # titles = check_reverse()
    t = len(titles)

    # 히든 매트릭스 차원: t(유전자 개수) * (n-k+1)
    # 프로필 매트릭스 차원: 4(염기) * k
    hidden_matrix = [[0 for _ in range((n-k+1)*2)] for _ in range(t)]
    profile_matrix = [[0 for _ in range(k)] for _ in range(4)]


    # 초기 프로필, 히든 만들기
    init_matrix(profile_matrix, hidden_matrix, titles, offset_hashmap, input_file, k, E)
    motifs = find_motifs(hidden_matrix, input_file, titles, offset_hashmap, k)
    print('initial motifs', motifs)
    # print("[INFO]: Creating initial motifs")
        
    idx = 0
    while True:
        profile_matrix = update_profile(hidden_matrix, profile_matrix, titles, offset_hashmap, input_file, k, E)

        new_motifs = find_motifs(hidden_matrix, input_file, titles, offset_hashmap, k)
        flag = True
        for i in range(len(new_motifs)):
            if new_motifs[i] != motifs[i]:
                flag = False
                break
        motifs = new_motifs

        # if flag:
        #     break
        hidden_matirx = update_hidden(hidden_matrix, profile_matrix, titles, offset_hashmap, input_file, k)
        idx += 1

        if idx > 5:
            break


    print(titles)
    print('[INFO]: EM algorithm is done.')
    print('[INFO]: Consensus motif is', motifs)
    print('[INFO]: iteration: ', idx)
