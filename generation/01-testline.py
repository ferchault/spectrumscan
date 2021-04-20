import itertools as it
import functools
import igraph
import sys
import tqdm
class IntegerPartitions(object):
    @staticmethod
    @functools.lru_cache(maxsize=64)
    def _do_partition(total, maxelements, around=None, maxdz=None):
        """ Builds all integer partitions of *total* split into *maxelements* parts.
		Note that ordering matters, i.e. (2, 1) and (1, 2) are district partitions. Moreover, elements of zero value are allowed. In all cases, the sum of all elements is equal to *total*.
		There is no guarantee as for the ordering of elements.
		If a center *around* is given, then a radius *maxdz* is required.
		Only those partitions are listed where the L1 norm of the distance between partition and *around* is less or equal to *maxdz*.
		Args:
			total:			The sum of all entries. [Integer]
			maxelements:	The number of elements to split into. [Integer]
			around:			Tuple of N entries. Center around which partitions are listed. [Integer]
			maxdz:			Maximum absolute difference in Z space from center *around*. [Integer]
		Returns:
			A list of all partitions as lists.
		"""
        if (around is None) != (maxdz is None):
            raise ValueError("Cannot define center or radius alone.")

        if maxelements == 1:
            if around is not None and maxdz < abs(total - around[-maxelements]):
                return []
            else:
                return [[total]]
        res = []

        # get range to cover
        if around is None:
            first = 0
            last = total
            limit = None
        else:
            first = max(0, around[-maxelements] - maxdz)
            last = min(total, around[-maxelements] + maxdz)
        for x in range(first, last + 1):
            if around is not None:
                limit = maxdz - abs(x - around[-maxelements])
            for p in IntegerPartitions._do_partition(
                total - x, maxelements - 1, around, limit
            ):
                res.append([x] + p)
        return res

    @staticmethod
    def partition(total, maxelements, around=None, maxdz=None):
        """ Builds all integer partitions of *total* split into *maxelements* parts.
		Note that ordering matters, i.e. (2, 1) and (1, 2) are district partitions. Moreover, elements of zero value are allowed. In all cases, the sum of all elements is equal to *total*.
		There is no guarantee as for the ordering of elements.
		If a center *around* is given, then a radius *maxdz* is required.
		Only those partitions are listed where the L1 norm of the distance between partition and *around* is less or equal to *maxdz*.
		Args:
			total:			The sum of all entries. [Integer]
			maxelements:	The number of elements to split into. [Integer]
			around:			Iterable of N entries. Center around which partitions are listed. [Integer]
			maxdz:			Maximum absolute difference in Z space from center *around*. [Integer]
		Returns:
			A list of all partitions as lists.
		"""
        if around is not None:
            return IntegerPartitions._do_partition(
                total, maxelements, tuple(around), maxdz
            )
        else:
            return IntegerPartitions._do_partition(total, maxelements)
def integer_partition(nsites, total):
    res = []
    for c in IntegerPartitions.partition(total, nsites):
        if max(c) > 2 :
            continue
        if c[0] == 0:
            continue
        if c[::-1] not in res:
            res.append(c)
    return res

def filtered_clean(bonded, partition, nhs, colors):
    return filtered_clean_cached(tuple(bonded), tuple(partition), nhs, colors)

@functools.lru_cache(maxsize=100)
def filtered_clean_cached(bonded, partition, nhs, colors):
    return filter_list(clean_list(bonded, partition, nhs), colors)

def clean_list(bonded, partition, nhs):
    res = []
    if len(partition) == 0:
        return [[]]
    for add_here_A in range(0, partition[0]+1):
        fill_A = partition[0] - add_here_A
        if fill_A > nhs:
            continue
        for added_A in it.combinations(bonded, add_here_A):
            group_A = list(added_A) + ['H'] * fill_A
            
            remain_bonded = set(bonded) - set(added_A)
            remain_partition = partition[1:]
            remain_nhs = nhs - fill_A
            for partial in clean_list(remain_bonded, remain_partition, remain_nhs):
                if len(partial) > 0 and (partial[0] == [] and group_A == []):
                    continue
                res.append([group_A] + partial)
    return res

def filter_list(clean, colors):
    accepted = []
    representations = []
    for setup in clean:
        representation = '-'.join(['.'.join(_) for _ in setup])
        for idx, color in enumerate(colors):
            if color == "1":
                representation = representation.replace(str(idx), 'O')
        if representation not in representations:
            accepted.append(setup)
            representations.append(representation)
    return accepted

def find_terminal_oxygens(edges, colors, nhs):
    terminal = []
    for aidx in range(len(colors)):
        if colors[aidx] == "0":
            continue
        nbonds = len([_ for _ in edges if _ == str(aidx)])
        if nbonds == 2:
            continue
        if nhs[aidx] == "1":
            continue
        terminal.append(str(aidx))
    return terminal

def verify_ring(decorations, other_ring_index, other_index_order, terminal):
    orders = []
    for site in range(6):
        order = 2
        for entry in decorations[site]:
            if entry == "H":
                order -= 1
                continue
            if entry == str(other_ring_index):
                order -= other_index_order
                continue
            if entry in terminal:
                order -= 2
                continue
            order -= 1
        orders.append(str(order))
    orders = ''.join(orders)
    #print (orders)
    return orders in ['000000', '110000', '011000', '001100', '000110', '000011', '100001', '111100', '011110', '001111', '100111', '110011', '111001', '111111', '121000', '012100', '001210', '000121', '100012','210001']

def permutate_mol(line):
    parts = line.replace("(", "").replace(")", "").replace(",", "").strip().split()[2:]
    colors, edges, nhs = tuple(parts[:9]), parts[9:-9], parts[-9:]
    terminal = find_terminal_oxygens(edges, colors, nhs)
    
    def get_bonded(ring):
        bonded = []
        for atoms in zip(edges[::2], edges[1::2]):
            if ring in atoms:
                bonded.append([_ for _ in atoms if _ != ring][0])
        return bonded
    
    ring1idx = colors.index("0")
    connect1 = get_bonded(str(ring1idx))
    ring2idx = len(colors)-1-colors[::-1].index("0")
    connect2 = get_bonded(str(ring2idx))
    
    njoins1 = len(connect1) + int(nhs[ring1idx])
    njoins2 = len(connect2) + int(nhs[ring2idx])
    
    elements = [8]*7 + [6] * 12 + [1]*12
    accepted = []
    dumps = []
    for partition1 in integer_partition(6, njoins1):
        counter = 0
        debug = 0
        thisaccepted = []

        print(len(clean_list(connect1, partition1, int(nhs[ring1idx]))))
        print(clean_list(connect1, partition1, int(nhs[ring1idx]))[:10])
        print(len(filter_list(clean_list(connect1, partition1, int(nhs[ring1idx])), colors)))

        return []
        for i in filter_list(clean_list(connect1, partition1, int(nhs[ring1idx])), colors):
            if not (verify_ring(i, ring2idx, 1, terminal) or verify_ring(i, ring2idx, 2, terminal)):
                continue

            for partition2 in integer_partition(6, njoins2):
                for j in filtered_clean(connect2, partition2, int(nhs[ring2idx]), colors):
                    test1 = (verify_ring(i, ring2idx, 1, terminal) and verify_ring(j, ring1idx, 1, terminal))
                    test2 = (verify_ring(i, ring2idx, 2, terminal) and verify_ring(j, ring1idx, 2, terminal))
                    if not (test1 or test2):
                        continue

                    #this_graph = get_graph(i, j, edges, ring1idx, ring2idx, nhs)
                    #for a in thisaccepted:
                    #    if this_graph.isomorphic_vf2(a, color1=elements, color2=elements):
                    #        break
                    #else:
                    #    thisaccepted.append(this_graph)
                    counter += 1

        print (counter)
        
        accepted.append(thisaccepted)
    return accepted
def get_graph(i,j, edges, ring1idx, ring2idx, nhs):
    bonded1, bonded2 = i, j
    
    oxygens = list(range(7))
    carbons = list(range(7,7+12))
    hydrogens = list(range(7+12, 7+12+12))
    
    # mapping
    mapping = dict()
    nrings = 0
    for nauty in range(9):
        if nauty in (ring1idx, ring2idx):
            mapping[str(nauty)] = carbons[nrings*6:(nrings+1)*6]
            nrings += 1
        else:
            mapping[str(nauty)] = [oxygens.pop(0)]    
    # BONDS
    bonds = []
    
    # rings
    for i, j in ((0,1), (1,2),(2,3), (3,4), (4,5), (5,0)):
        bonds.append(((i+7), (j+7)))
        bonds.append(((i+7+6), (j+7+6)))

    # atom-atom bonds
    for a, b in zip(edges[::2], edges[1::2]):
        a = mapping[a]
        b = mapping[b]
        
        if len(a) > 1 or len(b) > 1:
            continue
        bonds.append((a[0], b[0]))
        
    # OH bonds
    for idx, nh in enumerate(nhs):
        if idx not in (ring1idx, ring2idx) and nh != "0":
            bonds.append((mapping[str(idx)][0], hydrogens.pop(0)))
    
    # ring attached
    offset = 7
    ringring = []
    for ringgroups in (bonded1, bonded2):
        for site, toadd in enumerate(ringgroups):
            for entry in toadd:
                if entry == "H":
                    other = hydrogens.pop(0)
                else:
                    other = mapping[entry]
                    if len(other) == 1:
                        other = other[0]  
                    else:
                        ringring.append(offset+site)
                        continue
                bonds.append((offset+site, other))
        offset += 6
    
    assert(len(hydrogens) == 0)
    
    # ring-ring
    if len(ringring) == 2:
        bonds.append(tuple(ringring))
    
    graph = igraph.Graph(bonds)
    assert(len(graph.components()))
    #return bonds
    return graph

for line in open(sys.argv[1]):
    print ('#IN', line.strip())
    accepted = permutate_mol(line)
    for graph in accepted:
        print ('#OUT', ' '.join(['-'.join([str(__) for __ in _.tuple]) for _ in graph.es]))
