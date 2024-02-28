import time
from collections import Counter
from functools import reduce
from math import comb, factorial
from operator import mul

n = 17
max_arcs = [comb(i, 2) for i in range(0, n + 1)]
disconnected = [
    0,
    0,
    1,
    4,
    26,
    296,
    6064,
    230896,
    16886864,
    2423185664,
    687883494016,
    387139470010624,
    432380088071584256,
    959252253993204724736,
    4231267540316814507357184,
    37138269572860613284747227136,
    649037449132671895468073434916864,
    22596879313063804832510513481261154304,
    1568021213156739677613469882713171817332736,
    216941651871139111112257250053688796883675316224,
    59863107263844702400145424148028259703877851204222976,
    32954784141164432483664757174340608734349306273575787298816,
    36201100270287457206439712419041081799874552386281234238766055424,
    79370082993650493456763890642792332552505079606953465748596793392234496,
    347376267718609988990890055784978112019812340642745454746390010290681465012224,
    3035420144118576189529829298097903070436795979406208784555665386583711850333579771904,
    52962935384771904339216029368932375592288909108214621340719992253309601740739355993792577536,
    1845492799036376605599623847952100832765704496104863816834539882791539449497887801190970769385979904,
    128435922495457351369587905914421958994903392241439658360267117007336100805942319214416009856031712750862336,
    17854034057741747907609379605109267709693965827154209940879037366855540353557358894352884617265866096880888039079936,
    4957919765926901738601844851510927154063991826697182580924993402086705655705682061559710603812518801409345943179066588790784,
    2750488336567853875455808693604521769473058416369018727081258299239258257457843982132781555411441256143316922391990938000265341566976,
    3048582568667962802877924981092473001359916945835225583678419957521059771801285210753493386244882312181903820233608964302802213767524326047744,
    6751368128785793575074260524130484547319422123450346133227736696838823226140946255191458471790959774983471425936046950541918273604759428359476282916864,
    29875599416888414619771211712977950784970060838073223069962620168890598736375623231285529333676963821675502045627372059306654933089020431047407421088751593455616,
    264177369737507897334951631472063030259487672861334504106034435540784489608731246059241241996550516644727950351124531714015781186144443418233464101628917059290077922852864,
    4668205014991117090874210890974320131459807526502759206929993706935856097709373439774914694000315760392709274590442985276617489875206538767806512044916333114159836298348690392743936
]
connected = [1, 1] + [pow(2, comb(i, 2)) - disconnected[i] for i in range(2, len(disconnected))]
known = {
    (3, 2): 3,
    (3, 3): 1
}


# This is the structure for C(n,k)
class SubGraph:
    def __init__(self, size: int, arcs: int):
        self.size: int = size
        self.arcs: int = arcs

    def __str__(self):
        return f"C({self.size}, {self.arcs})"

    # Calculate C(n,k)
    def calculate(self) -> int:
        if self.size <= 2:
            return 1
        elif (self.size, self.arcs) in known.keys():
            return known[(self.size, self.arcs)]
        # print(str(self))

        result = 0
        base = comb(max_arcs[self.size], self.arcs)
        sequences = find_sequences(self.size, self.arcs)
        for seq in sequences:
            configs = find_configs(
                len(seq),
                self.arcs,
                self.arcs,
                [x - 1 for x in seq],
                [max_arcs[x] for x in seq]
            )
            patterns = [Pattern([SubGraph(seq[i], config[i]) for i in range(0, len(seq))]) for config in
                        configs]
            result += sum([x.calculate() for x in patterns]) * find_combinations(seq)
        known[(self.size, self.arcs)] = base - result
        return known[(self.size, self.arcs)]


# This is a list of SubGraph instances
class Pattern:
    def __init__(self, subgraphs: list[SubGraph]):
        self.subgraphs: list[SubGraph] = subgraphs

    def __str__(self):
        return ' * '.join([x.__str__() for x in self.subgraphs])

    # Calculating productivity of C(n,k)
    def calculate(self) -> int:
        return reduce(mul, [x.calculate() for x in self.subgraphs], 1)


# Calculate H(X1,...,Xk)
def find_combinations(sequence) -> int:
    combinations = 1
    start = sum(sequence)
    c = 0
    while start > 1:
        combinations *= comb(start, sequence[c])
        start -= sequence[c]
        c += 1
    duplicates = {key: value for key, value in Counter(sequence).items() if value > 1}
    for duplicate in duplicates.values():
        combinations //= factorial(duplicate)
    return combinations


# Calculate the set S
def find_sequences(n, min_arc=None, upper=None, seq=None, sequences=None) -> list[list[int]]:
    if seq is None:
        upper = n
        seq = []
        sequences = []
    if min_arc is None:
        min_arc = n - 1

    if n == 0:
        arcs_sum = sum([max_arcs[i] for i in seq])
        if len(seq) > 1 and arcs_sum >= min_arc:
            sequences.append(seq)
    else:
        for num in reversed(range(1, min(n, upper) + 1)):
            find_sequences(n - num, min_arc, num, seq + [num], sequences)
    return sequences


# Calculate the set T
def find_configs(n, k1, k2, lo_limits, hi_limits) -> list[list[int]]:
    # params = '-'.join([str(n), str(k1), str(k2), ''.join(str(lo_limits)), ''.join(str(hi_limits))])
    # if params in find_configs_cache.keys():
    #     return find_configs_cache[params]

    new_k2 = 9999 if k2 <= -1 else k2
    current_list = []

    if n == 1:
        return [[i] for i in range(max(lo_limits[0], k1), min(hi_limits[0], new_k2) + 1)]
    tails = find_configs(
        n - 1,
        k1 - hi_limits[0],
        k2 - lo_limits[0],
        lo_limits[1:],
        hi_limits[1:]
    )
    for tail in tails:
        for i in range(max(lo_limits[0], k1 - sum(tail)),  # max(3, 5) = 5
                       min(hi_limits[0], new_k2 - sum(tail)) + 1):  # min(6, inf) = 6
            current_list.append([i] + tail)

    # find_configs_cache[params] = current_list
    return current_list


# Calculate (X1,...,Xk)k1,k2
def calculate_seq(seq, disable_tricks=False):
    if seq == [n - 1, 1] and not disable_tricks:
        print("idea1")
        return (connected[n - 1] - SubGraph(n - 1, n - 2).calculate()) * n
    elif seq == [n - 2, 2] and not disable_tricks:
        print("idea2")
        return (connected[n - 2] - SubGraph(n - 2, n - 3).calculate()) * comb(n, n - 2)

    configs = find_configs(
        len(seq),
        n - 1,
        -1,
        [x - 1 for x in seq],
        [max_arcs[x] for x in seq]
    )
    patterns = [Pattern([SubGraph(seq[i], config[i]) for i in range(0, len(seq))]) for config in configs]
    results = [pattern.calculate() for pattern in patterns]
    return sum(results) * find_combinations(seq)


if __name__ == '__main__':
    start_time = time.time()
    sequences = find_sequences(n)
    disconnected = sum([comb(comb(n, 2), k) for k in range(0, n - 1)])
    counter = 0
    for seq in sequences:
        res = calculate_seq(seq)
        disconnected += res
        # print(seq, res)
        counter += 1
        print(counter * 10000 // len(sequences) / 100, '%')

    print("n =", n, "- disconnessi-->", disconnected)
    print("n =", n, "- connessi-->", pow(2, comb(n, 2)) - disconnected)
    if len(connected) == n:
        connected.append(pow(2, comb(n, 2)) - disconnected)
    print(f"Elapsed time: {(time.time() - start_time)} sec")
