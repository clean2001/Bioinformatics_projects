import sys, time
import enum
import random as rd

E = 0.02 # 엡실론

class Nucleotide(enum.Enum):
    A = 0
    C = 1
    G = 2
    T = 3

def open_and_parse(input_file: str):
    f = open(input_file, 'r')

    offset_hashmap, titles = dict(), list()
    max_len, seq, seqs= 0, '', list()
    while True:
        line = f.readline().rstrip()

        if not line:
            if len(seq):
                seqs.append(seq)

            f.close()
            break
        
        if line[0] == '>':
            length = len(seq)
            max_len = max(length, max_len)
            if len(seq):
                seqs.append(seq)
            seq = ''

            title = line[1:]
            offset_hashmap[title] = f.tell()
            titles.append(title)
        else:
            seq += line
        
    return titles, offset_hashmap, max_len, seqs



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


def update_hidden(seqs: list, rev_seqs:list, hidden: list, profile: list, motif_len: int):
    for i in range(len(seqs)):
        seq = seqs[i]
        for j in range(len(seq)-motif_len+1):
            probability = 1
            for k in range(motif_len):
                if seq[k+j] == 'A':
                    probability *= profile[Nucleotide.A.value][k]
                elif seq[k+j] == 'C':
                    probability *= profile[Nucleotide.C.value][k]
                elif seq[k+j] == 'G':
                    probability *= profile[Nucleotide.G.value][k]
                elif seq[k+j] == 'T':
                    probability *= profile[Nucleotide.T.value][k]
                else:
                    print('[Error]')
                    exit(-1)
            hidden[i][j] = probability

    for i in range(len(rev_seqs)):
        seq = rev_seqs[i]
        for j in range(len(seq)-motif_len+1):
            probability = 1
            for k in range(motif_len):
                if seq[k+j] == 'A':
                    probability *= profile[Nucleotide.A.value][k]
                elif seq[k+j] == 'C':
                    probability *= profile[Nucleotide.C.value][k]
                elif seq[k+j] == 'G':
                    probability *= profile[Nucleotide.G.value][k]
                elif seq[k+j] == 'T':
                    probability *= profile[Nucleotide.T.value][k]
                else:
                    print('[Error]')
                    exit(-1)
            hidden[i][j+len(seq)-motif_len+1] = probability
        
        total = sum(hidden[i])
        for j in range(2*(len(seq)-motif_len+1)):
            hidden[i][j] /= total
    
    return hidden


def update_profile(seqs, rev_seqs, hidden: list, profile: list, motif_len: int, eps:float):
    # 프로필에 엡실론 값 더해주기
    for row in profile:
        for j in range(len(row)):
            row[j] = eps

    for i in range(len(seqs)):
        seq = seqs[i]
        for j in range(len(seq)-motif_len+1):
            for k in range(motif_len):
                if seq[j+k] == 'A':
                    profile[Nucleotide.A.value][k] += hidden[i][j]
                elif seq[j+k] == 'C':
                    profile[Nucleotide.C.value][k] += hidden[i][j]
                elif seq[j+k] == 'G':
                    profile[Nucleotide.G.value][k] += hidden[i][j]
                elif seq[j+k] == 'T':
                    profile[Nucleotide.T.value][k] += hidden[i][j]
                else:
                    print('[Error]')
                    exit(-1)
    for i in range(len(rev_seqs)):
        seq = rev_seqs[i]
        col_num = len(seq) - motif_len + 1

        for j in range(col_num):
            for k in range(motif_len):
                if seq[j+k] == 'A':
                    profile[Nucleotide.A.value][k] += hidden[i][j+col_num]
                elif seq[j+k] == 'C':
                    profile[Nucleotide.C.value][k] += hidden[i][j+col_num]
                elif seq[j+k] == 'G':
                    profile[Nucleotide.G.value][k] += hidden[i][j+col_num]
                elif seq[j+k] == 'T':
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
    return profile

# 현재 프로필과 과거 프로필이 같은지를 확인
def is_same_profile(cur_profile, prev_profile):
    for i in range(len(cur_profile)):
        for j in range(len(cur_profile[i])):
            if cur_profile[i][j] != prev_profile[i][j]:
                return False
    
    return True


def init_matrix(profile:list, hidden:list, motif_len, eps: float, seqs: list, rev_seqs: list):
    # 프로필에 엡실론 값 더해주기
    for row in profile:
        for j in range(len(row)):
            row[j] = eps

    for i in range(4): # 랜덤한 수를 배정
        for j in range(motif_len):
            profile[i][j] = rd.random()
        
    for col in range(motif_len):
        total = 0
        
        for row in range(4):
            total += profile[row][col]

        for i in range(4):
            profile[i][col] /= total
    
    # hidden matrix 만들기
    hidden = update_hidden(seqs, rev_seqs, hidden, profile, motif_len)

    return profile, hidden


def find_motifs(seqs, rev_seqs, hidden, motif_length):
    motifs_array = []
    for i in range(len(seqs)):
        seq = seqs[i]
        col_num = len(seq) - motif_length + 1
        max_idx, max_prob = 0, 0
        for j in range(2*col_num):
            if max_prob < hidden[i][j]:
                max_idx = j
                max_prob = hidden[i][j]
        
        motif = ''
        if max_idx >= col_num:
            motif = rev_seqs[i][max_idx-col_num: max_idx+motif_length-col_num]
        else:
            motif = seqs[i][max_idx: max_idx+motif_length]

        motifs_array.append(motif)


    return motifs_array


if __name__=='__main__':
    start_time = time.time()
    if len(sys.argv) < 3:
        print("[ERROR]: Insufficient args")
        sys.exit()
    
    input_file = sys.argv[1]
    k = int(sys.argv[2]) # 보통 11, 12 mer
    titles, offset_hashmap, n, seqs = open_and_parse(input_file) # n은 max len
    rev_seqs = list()
    for seq in seqs:
        rev_seqs.append(process_reverse(seq))
    
    t = len(titles)


    # 히든 매트릭스 차원: t(유전자 개수) * (n-k+1)
    # 프로필 매트릭스 차원: 4(염기) * k
    hidden_matrix = [[0 for _ in range(2*(n-k+1))] for _ in range(t)]
    profile_matrix = [[0 for _ in range(k)] for _ in range(4)]


    # 초기 프로필, 히든 만들기
    init_matrix(profile_matrix, hidden_matrix, k, E, seqs, rev_seqs)
    motifs = find_motifs(seqs, rev_seqs, hidden_matrix, k)
    print('initial motifs', motifs)

        
    idx = 0
    while True:
        profile_matrix = update_profile(seqs, rev_seqs, hidden_matrix, profile_matrix, k, E)
        hidden_matirx = update_hidden(seqs, rev_seqs, hidden_matrix, profile_matrix, k)
        new_motifs = find_motifs(seqs, rev_seqs, hidden_matrix, k)
        flag = True
        for i in range(len(new_motifs)):
            if new_motifs[i] != motifs[i]:
                flag = False
                break
        if flag:
            break

        motifs = new_motifs

        idx += 1

        print("[INFO]: iteration ", idx)
        if idx > 10:
            break

    print('[INFO]: EM algorithm is done.')
    print('[INFO]: Final motifs are')
    for m in motifs:
        print(m)
    print('[INFO]: Total iteration: ', idx)

    end_time = time.time()

    print('[INFO]: Total execution time: ', end_time - start_time)