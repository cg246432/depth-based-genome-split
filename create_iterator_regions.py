
import glob
import collections

def coverage_across_chromosomes(cov_dir):
    global_mean_depth = []
    global_mean_depth_region_count = 0
    out = collections.defaultdict(dict)
    for i in glob.glob(f"{cov_dir}/*.cov"):
        chrom = i.replace(".cov", "").rsplit("/", 1)[1]
        with open(i, 'r') as f:
            header = f.readline().strip().split("\t")
            for row in f:
                row_dict = dict(zip(header, row.strip().split("\t")))
                out[chrom][row_dict['startpos']] = float(row_dict['meandepth'])
                global_mean_depth.append(float(row_dict['meandepth']))
                global_mean_depth_region_count += 1
    return out, (sum(global_mean_depth) / global_mean_depth_region_count)

def target_mean_depth(meandepth_dict):
    # Per chromosome mean depth
    return sum([v for v in list(meandepth_dict.values())]) / len(meandepth_dict)

def sum_of_differences(group_means, target_mean):
    # Value we're optimizing for
    # return sum([(i - target_mean) for i in group_means])
    return (sum([i for i in group_means])) - target_mean

def merge_indices(depths, indices):
    # Return a list of len n -1 from the input list where the two indices found are merged
    out_depths = []
    # print(indices)
    merged = depths[indices[0]] + depths[indices[1]]
    for i in range(0, len(depths)):
        if i == indices[0]:
            out_depths.append(merged)
        elif i == indices[1]:
            pass
        else:
            out_depths.append(depths[i])
    return out_depths

def search_for_smallest_diff(depths, target_mean_depth): # List of lists
    # Want to look for smallest/most negative contiguous merge relative to the target mean depth
    min_indices = [] # List with len 2, captures group indices to merge
    min_diff = None

    for i in range(1, len(depths)):
        d = sum_of_differences(depths[i] + depths[i-1], target_mean_depth)
        if not min_diff or d < min_diff:
            min_indices = [i-1, i]
            min_diff = d

    return merge_indices(depths, min_indices), min_indices

        
def merge_regions_by_differences(cov_dict, chrom, target_num_groups):
    chrom_mean = target_mean_depth(cov_dict[chrom])
    regions_in_order = [k for k in list(cov_dict[chrom].keys())]
    init_depths = [[v] for v in list(cov_dict[chrom].values())] # List of lists with size 1
    
    current_depths = init_depths
    while len(current_depths) > target_num_groups:
        current_depths, min_indices = search_for_smallest_diff(current_depths, chrom_mean)

    return current_depths, chrom_mean

def convert_final_depths_to_regions(cov_dict, chrom, target_num_groups):
    regions_in_order = collections.deque([k for k in list(cov_dict[chrom].keys())])
    chrom_depths, chrom_mean = merge_regions_by_differences(cov_dict, chrom, target_num_groups)
    out_regions_list = []
    for depth_sublist in chrom_depths:
        region_sublist = []
        len_region = len(depth_sublist)
        # While regions remain in the subgroup, return the exact number -- ultimately it won't be used except the first
        while len_region > 0:
            region_sublist.append(regions_in_order.popleft())
            len_region = len_region - 1
        # If there is exactly 1 region left, we want that region to be 1 BP less than the next region -- need to try looking ahead
        if len_region == 0:
            # Get the next number, the first in regions_in_order (first of the next subregion) and subtract 1 bp to make contiguous
            try:
                region_sublist.append(str(int(regions_in_order[0]) - 1))
            # This doesn't work if it's the last region of the chromosome
            except IndexError:
                # Want to let the script know that this is the end, and this number should ultimately be left off
                region_sublist.append("END")

        out_regions_list.append(region_sublist)
    return out_regions_list


def find_largest_diff(chrom_depth_dict):
    largest_val = None
    key = ""
    for k, v in chrom_depth_dict.items():
        if not largest_val or v > largest_val:
            key = k
            largest_val = v
    return key


def allocate_groups_per_chromosome(chrom_depth_dict, sample_wide_mean, total_groups=300):
    regions_to_allocate = total_groups
    regions_per_chrom = {}
    # Initialize the count
    for k in chrom_depth_dict.keys():
        regions_per_chrom[k] = 1
        regions_to_allocate = regions_to_allocate - 1

    current_chrom_depth_dict = chrom_depth_dict 
    while regions_to_allocate != 0:
        chrom = find_largest_diff(current_chrom_depth_dict)
        regions_per_chrom[chrom] += 1
        current_chrom_depth_dict[chrom] -= sample_wide_mean
        regions_to_allocate = regions_to_allocate - 1

    return regions_per_chrom

def format_regions(grouped_regions):
    out = []
    for subregion in grouped_regions:
        if len(subregion) > 1:
            out_string = f"{subregion[0]}-{subregion[-1]}"
        else:
            out_string = subregion[0]
        out_string = out_string.replace("-END", "") # samtools/htslib will interpret as start -> end of chromosome
        out.append(out_string)
    return out


def main():
    cov_dict, sample_wide_mean_depth = coverage_across_chromosomes("cov")
    # ROI = parse_roi()

    chrom_depth_dict = {}
    for k, depthdict in cov_dict.items():
        chrom_depth_dict[k] = target_mean_depth(depthdict)
    
    regions_per_chrom = allocate_groups_per_chromosome(chrom_depth_dict, sample_wide_mean_depth)
    for k in cov_dict.keys():
        grouped_regions = convert_final_depths_to_regions(cov_dict, k, regions_per_chrom[k])
        print(format_regions(grouped_regions))

if __name__=="__main__":
    main()
