import sys

def get_mid_distance(output, bed, distance = 1000, up_bin_count=25, down_bin_count=25,tx_bin_count=100, ):
    from collections import OrderedDict
    from statistics import mean
    total_number_of_bins = tx_bin_count + up_bin_count * 2
    pre_chr = ""
    pre_key = ""
    pre_end = 0

    out = open(output, "w")
    dis_dict = OrderedDict()

    f=open(bed, "r")
    for line in f.readlines():
        [chrom, start, end, name, score, strand] = line.strip().split("\t")
        start = int(float(start))
        end = int(float(end))
        key = f"{chrom}:{start}-{end}"
        dis_dict[key] = [chrom, start, end, name, str(strand)]

        if chrom == pre_chr:
            dis = int(start-pre_end)/2
            dis_dict[key].append(dis)
            dis_dict[pre_key].append(dis)
        else:
            dis_dict[key].append(0)
        pre_chr = chrom
        pre_key = key
        pre_end = end
        
    for key, val in dis_dict.items():
        if len(val) < 7:
            val.append(0)

        chrom = val[0]
        tx_start = int(val[1])
        tx_end = int(val[2])
        tx_name = val[3]
        strand = val[4]
        tx_length = tx_end - tx_start

        up_length = val[-2]
        if up_length == 0:
            up_length = distance
            #up_bin_count = 0

        down_length = val[-1]
        if down_length == 0:
            down_length = distance
            #down_bin_count = 0                

        # upstream
        up_starts = []
        up_ends = []
        if up_bin_count > 0:
            up_bin_length = int(up_length / up_bin_count)
            up_start = tx_start - up_length
            up_starts.append(up_start)
            for i in range(up_bin_count - 1):
                up_starts.append(up_starts[-1] + up_bin_length)
            up_ends = up_starts[1:]
            up_ends.append(tx_start)

        # transcript
        ideal_bin_size = tx_length / tx_bin_count
        bin_size = int(ideal_bin_size)
        tx_starts = [tx_start]
        bin_sizes = [ideal_bin_size]
        for i in range(1, tx_bin_count):
            if mean(bin_sizes) < ideal_bin_size:
                tx_starts.append(tx_starts[-1] + int(bin_size) + 1)
                bin_sizes.append(tx_starts[-1] - tx_starts[-2])
            else:
                tx_starts.append(tx_starts[-1] + int(bin_size))
                bin_sizes.append(tx_starts[-1] - tx_starts[-2])
        tx_ends = tx_starts[1:]
        tx_ends.append(tx_end)

        # downstream
        down_starts = []
        down_ends = []
        if down_bin_count > 0:
            down_bin_length = int(down_length / down_bin_count)
            down_starts = [tx_end + x * down_bin_length for x in range(down_bin_count)]
            down_ends = [tx_end + (x + 1) * down_bin_length for x in range(down_bin_count)]

        starts = up_starts + tx_starts + down_starts
        ends = up_ends + tx_ends + down_ends

        total_number_of_bins = up_bin_count + down_bin_count + tx_bin_count
        for i in range(total_number_of_bins):
            out.write('{chrom}\t{tx_start}\t{tx_end}\t{tx_name}\t{bin_no}\t{strand}\n'.format(
                chrom=chrom,
                tx_start=int(starts[i]),
                tx_end=int(ends[i]),
                bin_no=(i + 1),
                tx_name=tx_name,
                strand=strand,
            ))

    out.close()

output = sys.argv[1]
bed = sys.argv[2]
distance = int(sys.argv[3])
get_mid_distance(output=output,bed=bed,distance=distance)

