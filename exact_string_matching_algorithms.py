import sys
import bisect
from collections import defaultdict

from io import StringIO
from math import ceil, log2


###########################################################################################
# Global variables
pattern1 = 'ATGATG'
pattern2 = 'CTCTCTA'
pattern3 = 'TCACTACTCTCA'

# Canis lupus familiaris genome, chromosome 1
file_name1 = 'cfa_ref_CanFam3.1_chr1.fa'
# Phoenix dactylifera genome
file_name2 = '42345_ref_DPV01_chrUn.fa'
# Ananascomosus genome, chromosome 1
file_name3 = '4615_ref_ASM154086v1_chr1.fa'


###########################################################################################
# Sorted index
class IndexSorted(object):

    # Constructor
    def __init__(self, name, text, length):
        self.name = name  # short sequence name
        self.text = text  # text to parse
        self.length = length  # length of substrings - index entries
        self.index = []  # list index

        for i in range(len(text) - length + 1):
            self.index.append((text[i:i + length], i))  # add <substr, offset> pair
        self.index.sort()  # sort pairs

    # Find possible pattern locations in index
    def query(self, p):
        st = bisect.bisect_left(self.index, (p[:self.length], -1))  # binary search
        en = bisect.bisect_right(self.index, (p[:self.length], sys.maxsize))  # binary search
        hits = self.index[st:en]  # this range of elements corresponds to the hits
        return [h[1] for h in hits]  # return just the offsets


###########################################################################################
# Hash table
class IndexHash(object):

    # Constructor
    def __init__(self, name, text, length):
        self.name = name  # short sequence name
        self.text = text  # text to parse
        self.length = length  # length of substrings -index entries
        self.index = {}  # dictionary index

        # Populate the dictionary with substring : [location1, location2, ...] pairs
        for i in range(len(text) - length + 1):
            substr = text[i:i + length]
            if substr in self.index:
                self.index[substr].append(i)  # substring already in dictionary
            else:
                self.index[substr] = [i]  # add to dictionary

    # Find possible pattern locations in dictionary
    def query(self, pattern):
        return self.index.get(pattern[:self.length], [])


###########################################################################################
# Suffix array
class SuffixArray(object):

    def sort_bucket(self, text, bucket, order):
        d = defaultdict(list)
        for i in bucket:
            key = text[i:i + order]
            d[key].append(i)
        result = []
        for k, v in sorted(d.items()):
            if len(v) > 1:
                result += self.sort_bucket(text, v, order * 2)
            else:
                result.append(v[0])
        return result

    def suffix_array(self):
        return self.sort_bucket(self.text, (i for i in range(len(self.text))), 1)

    def __init__(self, name, t):
        self.name = name  # short sequence name
        self.text = t
        self.index = self.suffix_array()

    def query(self, p):
        first = 0
        list = []
        last = len(self.index) - 1
        while first <= last:  # binary search
            midpoint = (first + last) // 2
            startIndex = self.index[midpoint]
            #print(self.text[startIndex:startIndex+len(p)])
            #print("startIndex: " + str(startIndex) + ", len(p): " + str(len(p)) + " , len(self.text): " + str(len(self.text)))
            if (startIndex + len(p)) <= len(self.text) and p == self.text[startIndex:startIndex + len(
                    p)]:  # search for matching first length(p) characters
                #print("match midpoint: " + str(midpoint))
                list.append(self.index[midpoint])
                j = midpoint - 1
                while True:  # find all matching before found string
                    if j < 0:
                        break
                    if (self.index[j] + len(p)) <= len(self.text) and p == self.text[self.index[j]:self.index[j] + len(p)]:
                        #print("match j: " + str(j))
                        list.append(self.index[j])
                        j = j - 1
                    else:
                        #print("break j: " + str(j))
                        break
                k = midpoint + 1
                while True:  # find all matching after found string
                    if k >= len(self.index):
                        break
                    if (self.index[k] + len(p)) <= len(self.text) and p == self.text[self.index[k]:self.index[k] + len(p)]:
                        #print("match k: " + str(k))
                        list.append(self.index[k])
                        k = k + 1
                    else:
                        #print("break k: " + str(k))
                        break
                return list
            else:
                if p < self.text[startIndex:startIndex + len(p)]:
                    #print("no match midpoint: " + str(midpoint) + ", first: " + str(first) + ", last: " + str(last))
                    last = midpoint - 1
                else:
                    #print("no match midpoint: " + str(midpoint) + ", first: " + str(first) + ", last: " + str(last))
                    first = midpoint + 1
        return list

###########################################################################################
# Suffix Tree
class _SNode():
    """Class representing a Node in the Suffix tree."""

    def __init__(self, idx=-1, parentNode=None, depth=-1):
        # Links
        self._suffix_link = None
        self.transition_links = []
        # Properties
        self.idx = idx
        self.depth = depth
        self.parent = parentNode
        self.generalized_idxs = {}

    def __str__(self):
        return ("SNode: idx:" + str(self.idx) + " depth:" + str(self.depth) +
                " transitons:" + str(self.transition_links))

    def _add_suffix_link(self, snode):
        self._suffix_link = snode

    def _get_suffix_link(self):
        if self._suffix_link != None:
            return self._suffix_link
        else:
            return False

    def _get_transition_link(self, suffix):
        for node, _suffix in self.transition_links:
            if _suffix == '__@__' or suffix == _suffix:
                return node
        return False

    def _add_transition_link(self, snode, suffix=''):
        tl = self._get_transition_link(suffix)
        if tl:  # TODO: imporve this.
            self.transition_links.remove((tl, suffix))
        self.transition_links.append((snode, suffix))

    def _has_transition(self, suffix):
        for node, _suffix in self.transition_links:
            if _suffix == '__@__' or suffix == _suffix:
                return True
        return False

    def is_leaf(self):
        return self.transition_links == []

    def _traverse(self, f):
        for (node, _) in self.transition_links:
            node._traverse(f)
        f(self)

    def _get_leaves(self):
        if self.is_leaf():
            return [self]
        else:
            return [x for (n, _) in self.transition_links for x in n._get_leaves()]


class SuffixTree():
    """Class representing the suffix tree."""

    def __init__(self, name, input):
        self.text = input
        self.name = name
        self.root = _SNode()
        self.root.depth = 0
        self.root.idx = 0
        self.root.parent = self.root
        self.root._add_suffix_link(self.root)

        self.build(input)

    def build(self, x):
        x += '$'
        """Builds a Suffix tree.
        Builds a Suffix tree using McCreight O(n) algorithm."""
        self.word = x
        u = self.root
        d = 0
        for i in range(len(x)):
            while u.depth == d and u._has_transition(x[d + i]):
                u = u._get_transition_link(x[d + i])
                d = d + 1
                while d < u.depth and x[u.idx + d] == x[i + d]:
                    d = d + 1
            if d < u.depth:
                u = self._create_node(x, u, d)
            self._create_leaf(x, i, u, d)
            if not u._get_suffix_link():
                self._compute_slink(x, u)
            u = u._get_suffix_link()
            d = d - 1
            if d < 0:
                d = 0

    def _create_node(self, x, u, d):
        i = u.idx
        p = u.parent
        v = _SNode(idx=i, depth=d)
        v._add_transition_link(u, x[i + d])
        u.parent = v
        p._add_transition_link(v, x[i + p.depth])
        v.parent = p
        return v

    def _create_leaf(self, x, i, u, d):
        w = _SNode()
        w.idx = i
        w.depth = len(x) - i
        u._add_transition_link(w, x[i + d])
        w.parent = u
        return w

    def _compute_slink(self, x, u):
        d = u.depth
        v = u.parent._get_suffix_link()
        while v.depth < d - 1:
            v = v._get_transition_link(x[u.idx + v.depth + 1])
        if v.depth > d - 1:
            v = self._create_node(x, v, d - 1)
        u._add_suffix_link(v)

    def query(self, y):
        node = self.root
        while True:
            edge = self._edgeLabel(node, node.parent)
            if edge.startswith(y):
                break

            i = 0
            while (i < len(edge) and edge[i] == y[0]):
                y = y[1:]
                i += 1

            if i != 0:
                if i == len(edge) and y != '':
                    pass
                else:
                    return []

            node = node._get_transition_link(y[0])
            if not node:
                return []

        leaves = node._get_leaves()
        return [n.idx for n in leaves]

    def _edgeLabel(self, node, parent):
        """Helper method, returns the edge label between a node and it's parent"""
        return self.word[node.idx + parent.depth: node.idx + node.depth] 


###########################################################################################
# Utility function for parsing input FASTA files
# Returns dictionary with sequence_name : base_string pairs
def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            if (current_short_name != None):
                fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.items():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa


###########################################################################################
# Pass through all possible alignments and leave only real ones
def filter_real_alignments(pattern, possible_alignments, text):
    real_alignments = []
    for i in possible_alignments:
        #print("text[" + str(i) + ":" + str(i+len(pattern)) + "]")
        if pattern == text[i:i + len(pattern)]:
            real_alignments.append(i)
    return real_alignments


###########################################################################################
# Choose string matching algorithm
def choose_algorithm():
    global choice
    global structure_name
    
    print('Exact String Matching Algorithms\n')
    print('Choose algorithm:')
    print('1. Sorted Index')
    print('2. Hash Table')
    print('3. Suffix Array')
    print('4. Suffix Tree')

    while True:
        try:
            choice = int(input())
            if (choice < 1 or choice > 4):
                raise ValueError
            
            if (choice == 1):
                structure_name = 'Sorted Index'
            elif (choice == 2):
                structure_name = 'Hash Table'
            elif (choice == 3):
                structure_name = 'Suffix Array'
            elif (choice == 4):
                structure_name = 'Suffix Tree'

            print('\nYou chose the ' + structure_name + ' algorithm.')

            return choice
        except ValueError:
            print('***ERROR***\n', 'Enter a valid number between 1 and 4:')


###########################################################################################
# Open file, create sequence dictionaries
def prepare_file(file_name):
    global parsed_fasta
    global file

    print('\nPreparing file ' + file_name + '...')
    file = open(file_name, 'r')
    text = file.read()
    string_io = StringIO(text)
    parsed_fasta = parse_fasta(string_io)
    print('Prepared!')


###########################################################################################
# Create adequate data structure, based on the chosen algorithm
def init_structure(name, text, length):
    if (choice == 1):
        return IndexSorted(name, text, length)
    elif (choice == 2):
        return IndexHash(name, text, length)
    elif (choice == 3):
        return SuffixArray(name, text)
    elif (choice == 4):
        return SuffixTree(name, text)
    else:
        return None

###########################################################################################
# Do all the needed processing for one file and given pattern
def per_pattern_processing(pattern):
    pattern_align_count = 0

    print('\nPattern ' + pattern)

    # Query each sequence
    for index in indexes:
        # Query table for possible alignment positions
        possible_alignments = index.query(pattern)

        # Find real alignment positions
        real_alignments = filter_real_alignments(pattern, possible_alignments, index.text)

        #print('[sequence: ' + index.name + ', pattern: ' + pattern + '], number of alignment positions: ' + str(
        #    len(real_alignments)))
        pattern_align_count += len(real_alignments)

    print('\nPattern ' + pattern + ' alignments in all sequences: ' + str(pattern_align_count))


###########################################################################################
# Create one data structure for each sequence
def create_indexes():
    # Array of index structures for each sequence in file
    global indexes
    indexes = []
    for key, value in parsed_fasta.items():
        print('Creating index for sequence ' + key)
        indexes.append(init_structure(key, value, 5))


###########################################################################################
# Do all the needed processing for one given file and three different patterns
def file_processing(file_name):
    prepare_file(file_name)

    create_indexes()

    per_pattern_processing(pattern1)
    per_pattern_processing(pattern2)
    per_pattern_processing(pattern3)

    file.close()


###########################################################################################
# Main function
def main():

    choose_algorithm()

    file_processing(file_name1)

    file_processing(file_name2)

    file_processing(file_name3)