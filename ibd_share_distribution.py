#plot IBD sharing(total cM) for pair of unrelated individuals randomly sampled in a well-mixed population
import numpy as np
import matplotlib.pyplot as plt
import gzip
from collections import defaultdict

def readIBD(ibdFile):
    ibd_length = defaultdict(0)
    with gzip.open(ibdFile, 'rt') as f:
        line = f.readline()
        while line:
            ind1, hap1, ind2, hap2, chr, start_bp, end_bp, len_cM = line.strip().split('\t')
            ind1, ind2 = sorted(ind1, ind2)
            pair = (ind1, ind2)
            ibd_length[pair] += float(len_cM)
            line = f.readline()
    return ibd_length


def main():
    parser = argparse.ArgumentParser(description='plot IBD sharing(total cM) for pair of unrelated individuals randomly sampled in a well-mixed population')
    parser.add_argument('--ibd', action="store", dest="ibd", type=str, required=True, help='path to ibd segment file.')
    #parser.add_arugment('-n', action="store", dest="num_samples", type=int, required=True, help="number of diploid samples")
    args = parser.parse_args()

    lens = np.array([value for key, value in readIBD(args.ibd).items()])
    plt.hist(lens, weights=np.zeros_like(lens)+1/lens.size)
    plt.xlabel('total IBD Sharing Length per pair of diploid individuals')
    plt.save('ibdSharing.png', dpi=300)


