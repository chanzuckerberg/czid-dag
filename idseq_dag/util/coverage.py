import math
import idseq_dag.util.log as log

# Turn (1, 5) into (0, 5) and (5, 1) into (5, 0).
def decrement_lower_bound(bounds):
    (bound_one, bound_two) = bounds
    if bound_one < bound_two:
        return (bound_one - 1, bound_two)
    else:
        return (bound_one, bound_two - 1)


def align_range(bounds):
    (bound_one, bound_two) = bounds
    return (min(bound_one, bound_two), max(bound_one, bound_two))


# Transform a range from one scale to another.
# The 'first_range' should be within the interval [first_scale_start, first_scale_end]
# For example, transform_range((3, 5) 0, 10, 100, 200) would return (130, 150)
def transform_range(
    first_range, first_scale_start, first_scale_end, second_scale_start, second_scale_end
):
    def transform_value(value):
        position = (value - first_scale_start) / (first_scale_end - first_scale_start)
        return position * (second_scale_end - second_scale_start) + second_scale_start

    return [transform_value(value) for value in first_range]

# One-off formatter that adds a couple more sig-figs to small numbers so they don't show up as "0.0"
def format_number(number):
    if number < 0.1:
        return round(number, 3)
    return round(number, 1)

def format_percent(number):
    return round(number, 3)

def round_if_within_epsilon(num, epsilon=0.001):
    r_num = round(num)

    if abs(r_num - num) < epsilon:
        return r_num
    else:
        return num

def safe_floor(num):
    return math.floor(round_if_within_epsilon(num))

def safe_ceil(num):
    return math.ceil(round_if_within_epsilon(num))

# Calculate the total distance covered by a set of overlapping range endpoints.
# Each range endpoint is [bp, depth_change], where bp is the base pair where the endpoint is,
# and depth_change is 1 for range start and -1 for range end.
def calculate_covered_length(endpoints):
    total_covered_length = 0
    cur_depth = 0
    last_start_point = 0

    endpoints.sort(key=lambda endpoint: endpoint[0])

    for endpoint in endpoints:
        cur_depth += endpoint[1]

        if cur_depth < 0:
            raise ValueError("coverage depth of %s is invalid. Malformed endpoints" % cur_depth)

        # Covered length is starting. Set the last start point.
        if endpoint[1] == 1 and cur_depth == 1:
            last_start_point = endpoint[0]
        # Covered length is ending. Add the covered distance to the total.
        elif endpoint[1] == -1 and cur_depth == 0:
            total_covered_length += endpoint[0] - last_start_point

    if cur_depth != 0:
        raise ValueError("coverage depth is %s after traversel. 0 is expected" % cur_depth)

    return total_covered_length


def calculate_accession_coverage(accession_id, accession_data, contigs_map, reads_map, num_bins):
    bin_size = accession_data["total_length"] / num_bins
    coverage = [{
        "depth": 0,
        "endpoints": [],
        "num_reads": 0,
        "num_contigs": 0
    } for i in range(num_bins)]

    for contig_name in accession_data["contigs"]:
        contig_obj = contigs_map[contig_name]
        # Ignore contigs with accession mismatch
        if contig_obj["accession"] != accession_id:
            continue

        # The bins and coverage array are 0-indexed, but subject start/end and coverage start/end are 1-indexed.
        # We convert everything to 0-index here and stay in 0-index for the rest of the code.
        # NOTE: We decrement only the lower bound here basically so that we can treat the discrete integer indices as a continuous interval.
        # We convert back to integer indices when we calculate coverage_arr_start/_end.
        (subject_start, subject_end) = decrement_lower_bound((contig_obj["subject_start"], contig_obj["subject_end"]))
        (query_start, query_end) = decrement_lower_bound((contig_obj["query_start"], contig_obj["query_end"]))

        # Find all bins that this contig overlaps, and calculate average coverage for each bin separately.
        (bin_start, bin_end) = align_range((subject_start / bin_size, subject_end / bin_size))

        for i in range(safe_floor(bin_start), safe_ceil(bin_end)):
            # Get the section of the accession that corresponds to the current bin and overlaps with the contig.
            accession_range = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

            # Convert the accession range to a contig range by using the alignment values.
            contig_range = transform_range(accession_range, subject_start, subject_end, query_start, query_end)

            # The contig coverage array should be the same length as the contig length.
            # If not, convert to the appropriate range in the coverage array.
            if contig_obj["total_length"] == len(contig_obj["coverage"]):
                coverage_range = align_range((contig_range[0], contig_range[1]))
            else:
                coverage_range = transform_range(contig_range, 0, contig_obj["total_length"], 0, len(contig_obj["coverage"]))
                coverage_range = align_range((coverage_range[0], coverage_range[1]))

            (coverage_arr_start, coverage_arr_end) = (safe_floor(coverage_range[0]), safe_ceil(coverage_range[1]))

            # Get the average coverage for the section of the contig that overlaps with this bin.
            avg_coverage_for_coverage_range = sum(contig_obj["coverage"][coverage_arr_start: coverage_arr_end]) / (coverage_arr_end - coverage_arr_start)

            # Multiply by the proportion of the bin that the contig covers.
            avg_coverage_for_bin = avg_coverage_for_coverage_range * (abs(accession_range[1] - accession_range[0]) / bin_size)

            coverage[i]["depth"] += avg_coverage_for_bin
            coverage[i]["endpoints"].append([max(i * bin_size, accession_range[0]), 1])
            coverage[i]["endpoints"].append([min((i + 1) * bin_size, accession_range[1]), -1])
            coverage[i]["num_contigs"] += 1

    # The processing for reads is similar to contigs above, but the avg coverage on the read is simply 1.
    for read_name in accession_data["reads"]:
        read_obj = reads_map[read_name]
        # Ignore contigs with accession mismatch
        if read_obj["accession"] != accession_id:
            continue

        # The bins and coverage array are 0-indexed, but subject start/end and coverage start/end are 1-indexed.
        # We convert everything to 0-index here and stay in 0-index for the rest of the code.
        # NOTE: We decrement only the lower bound here basically so that we can treat the discrete integer indices as a continuous interval.
        # We convert back to integer indices when we calculate coverage_arr_start/_end.
        (subject_start, subject_end) = decrement_lower_bound((read_obj["subject_start"], read_obj["subject_end"]))

        # Find all bins that this contig overlaps, and calculate average coverage for each bin separately.
        (bin_start, bin_end) = align_range((subject_start / bin_size, subject_end / bin_size))

        for i in range(safe_floor(bin_start), safe_ceil(bin_end)):
            # Get the section of the accession that corresponds to the current bin and overlaps with the read.
            accession_range = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

            # The read coverage is 1. Multiply by the proportion of the bin that the read covers.
            avg_coverage_for_bin = (abs(accession_range[1] - accession_range[0]) / bin_size)

            coverage[i]["depth"] += avg_coverage_for_bin
            coverage[i]["endpoints"].append([max(i * bin_size, accession_range[0]), 1])
            coverage[i]["endpoints"].append([min((i + 1) * bin_size, accession_range[1]), -1])
            coverage[i]["num_reads"] += 1

    final_coverage = []


    # For each index, an array of numbers is generated.
    # The array is sparse. Only bins with nonzero coverage are included.
    for index in range(len(coverage)):
        coverage_obj = coverage[index]
        if coverage_obj["depth"] > 0:
            final_coverage.append([
                index, # bin index
                format_number(coverage_obj["depth"]), # average coverage depth
                format_percent(calculate_covered_length(coverage_obj["endpoints"]) / bin_size), # coverage breadth
                coverage_obj["num_contigs"], # number of contigs
                coverage_obj["num_reads"], # number of reads
            ])

    return (final_coverage, bin_size)


# Calculate various accession stats.
def calculate_accession_stats(accession_data, contigs_map, reads_map):
    max_aligned_length = 0
    coverage_sum = 0
    endpoints = []
    prop_total_mismatch = 0

    for contig_name in accession_data["contigs"]:
        contig_obj = contigs_map[contig_name]

        (accession_start, accession_end) = align_range(decrement_lower_bound((contig_obj["subject_start"], contig_obj["subject_end"])))
        (contig_start, contig_end) = align_range(decrement_lower_bound((contig_obj["query_start"], contig_obj["query_end"])))

        accession_alignment_length = accession_end - accession_start
        if accession_alignment_length > max_aligned_length:
            max_aligned_length = accession_alignment_length

        # Restrict to the part of the coverage that corresponds to the alignment.
        coverage_sum += sum(contig_obj["coverage"][contig_start: contig_end])
        # TODO(mark): Should we multiply??
        prop_total_mismatch += contig_obj["prop_mismatch"]

        endpoints.append([accession_start, 1])
        endpoints.append([accession_end, -1])

    for read_name in accession_data["reads"]:
        read_obj = reads_map[read_name]

        (accession_start, accession_end) = align_range(decrement_lower_bound((read_obj["subject_start"], read_obj["subject_end"])))

        read_length = accession_end - accession_start
        if read_length > max_aligned_length:
            max_aligned_length = read_length

        coverage_sum += read_length
        prop_total_mismatch += read_obj["prop_mismatch"]

        endpoints.append([accession_start, 1])
        endpoints.append([accession_end, -1])

    covered_length = calculate_covered_length(endpoints)

    return {
        "max_aligned_length": max_aligned_length,
        "coverage_depth": coverage_sum / accession_data["total_length"],
        "coverage_breadth": covered_length / accession_data["total_length"],
        # Sum up the total prop mismatch and divide by number of contigs and reads.
        "avg_prop_mismatch": prop_total_mismatch / (len(accession_data["contigs"]) + len(accession_data["reads"]))
    }


def get_hit_group_json(contig_objs, read_objs, bin_size):
    num_contigs = len(contig_objs)
    num_reads = len(read_objs)
    contig_r = sum([contig_obj["num_reads"] for contig_obj in contig_objs])

    contig_byteranges = [contig_obj["byterange"] for contig_obj in contig_objs]

    # contig_obj and read_obj are very similar, except contig_obj has "num_reads".
    # From here onwards, we treat them the same.
    hit_objs = contig_objs + read_objs
    num_hits = num_contigs + num_reads

    def avg_field(field):
        return sum(map(lambda hit_obj: hit_obj[field], hit_objs)) / num_hits

    endpoints = []

    for hit_obj in hit_objs:
        endpoints.append(hit_obj["subject_start"])
        endpoints.append(hit_obj["subject_end"])

    hit_group_start = min(endpoints)
    hit_group_end = max(endpoints)

    hit_group_midpoint = ((hit_group_start - 1) + hit_group_end) / 2
    hit_group_bin_index = safe_floor(hit_group_midpoint / bin_size)

    return [
        num_contigs,
        num_reads,
        # Total number of contig reads.
        contig_r,
        # Alignment range
        hit_group_start,
        hit_group_end,
        # Alignment length. Can be different from alignment range.
        format_number(avg_field("alignment_length")),
        # Percent identity
        format_percent(avg_field("percent_id") / 100),
        # Number of mismatches
        format_number(avg_field("num_mismatches")),
        # Number of gaps
        format_number(avg_field("num_gaps")),
        # Bin index of midpoint. Used to display read group if read group is too small.
        hit_group_bin_index,
        # Byteranges in the contigs.fasta file for each contig.
        contig_byteranges
    ]


# In order to display numerous tiny hits, we organize them into groups.
def generate_hit_group_json(accession_data, accession_id, contigs_map, reads_map, num_bins):
    individual_reads = []
    individual_contigs = []

    bin_size = accession_data["total_length"] / num_bins
    read_bins = [[] for i in range(num_bins)]
    contig_bins = [[] for i in range(num_bins)]

    def process_hit(hit_type, hit_name):
        hit_map = reads_map if hit_type == "read" else contigs_map

        if not hit_name in hit_map:
            log.write(f"Could not find {hit_type} in map: {hit_name}")
            return

        hit_obj = hit_map[hit_name]

        # hitsummary is more strict than reassigned, since we judge whether a "hit" occurred.
        # Sometimes reassigned will have a value, but hitsummary won't.
        if hit_obj["accession"] != accession_id:
            log.write(f"Mismatched accession for {hit_name}: {hit_obj['accession']} (reassigned) versus {accession_id} (hitsummary)")
            return

        (accession_start, accession_end) = align_range(decrement_lower_bound((hit_obj["subject_start"], hit_obj["subject_end"])))

        # If the hit is larger than the bin size, treat it as an individual hit.
        if accession_end - accession_start >= bin_size:
            if hit_type == "read":
                individual_reads.append(hit_name)
            else:
                individual_contigs.append(hit_name)
        # Otherwise, put the hit into a bin based on its midpoint
        else:
            hit_midpoint = (accession_end + accession_start) / 2
            hit_bin_index = safe_floor(hit_midpoint / bin_size)

            if hit_type == "read":
                read_bins[hit_bin_index].append(hit_name)
            else:
                contig_bins[hit_bin_index].append(hit_name)

    for read_name in accession_data["reads"]:
        process_hit("read", read_name)

    for contig_name in accession_data["contigs"]:
        process_hit("contig", contig_name)

    hit_groups = []
    for read_name in individual_reads:
        read_obj = reads_map[read_name]
        hit_groups.append(get_hit_group_json([], [read_obj], bin_size))

    for contig_name in individual_contigs:
        contig_obj = contigs_map[contig_name]
        hit_groups.append(get_hit_group_json([contig_obj], [], bin_size))

    # Aggregate the stats of the binned reads.
    for i in range(num_bins):
        reads = read_bins[i]
        contigs = contig_bins[i]

        if len(reads) + len(contigs) == 0:
            continue
        else:
            read_objs = list(map(lambda read_name: reads_map[read_name], reads))
            contig_objs = list(map(lambda contig_name: contigs_map[contig_name], contigs))
            hit_groups.append(get_hit_group_json(contig_objs, read_objs, bin_size))

    return hit_groups


def generate_coverage_viz_json(accession_data, contigs_map, reads_map, max_num_bins):
    coverage_viz_json = {}

    for accession_id, accession_obj in accession_data.items():
        # Number of bins to calculate coverage for.
        num_bins = min(max_num_bins, accession_obj["total_length"])

        total_length = accession_obj["total_length"]
        hit_groups = generate_hit_group_json(accession_obj, accession_id, contigs_map, reads_map, num_bins)

        (coverage, coverage_bin_size) = calculate_accession_coverage(
            accession_id, accession_obj, contigs_map, reads_map, num_bins
        )

        accession_stats = calculate_accession_stats(
            accession_obj,
            contigs_map,
            reads_map
        )

        coverage_viz_json[accession_id] = {
            "total_length": total_length,
            "name": accession_obj["name"],
            "hit_groups": hit_groups,
            "coverage": coverage,
            "coverage_bin_size": coverage_bin_size,
            "max_aligned_length": accession_stats["max_aligned_length"],
            "coverage_depth": format_number(accession_stats["coverage_depth"]),
            "coverage_breadth": format_percent(accession_stats["coverage_breadth"]),
            "avg_prop_mismatch": format_percent(accession_stats["avg_prop_mismatch"]),
        }

    return coverage_viz_json


# Mutates input parameters.
def remove_taxons_with_no_contigs(accession_data, taxons_to_accessions):
    taxons_to_remove = []

    for taxon, accessions in taxons_to_accessions.items():
        total_contigs = sum([len(accession_data[accession_id]["contigs"]) for accession_id in accessions])

        if total_contigs == 0:
            taxons_to_remove.append(taxon)

    # Remove all taxons with no mapped contigs.
    for taxon in taxons_to_remove:
        accessions = list(taxons_to_accessions[taxon])
        del taxons_to_accessions[taxon]
        for accession_id in accessions:
            del accession_data[accession_id]


def get_unassigned_reads_set(accession_data):
    unassigned_reads_set = set()
    for _, accession_obj in accession_data.items():
        for read in accession_obj["reads"]:
            unassigned_reads_set.add(read)

    return unassigned_reads_set


def generate_coverage_viz_summary_json(taxons_to_accessions, accession_data, coverage_viz_obj, taxons_to_total_accession_count):
    coverage_viz_summary_json = {}

    for taxon, accessions in taxons_to_accessions.items():
        coverage_viz_summary_json[taxon] = {
            "best_accessions": list(map(lambda accession_id: {
                "id": accession_id,
                "name": accession_data[accession_id]["name"],
                "num_contigs": len(accession_data[accession_id]["contigs"]),
                "num_reads": len(accession_data[accession_id]["reads"]),
                "score": format_number(accession_data[accession_id]["score"]),
                "coverage_depth": coverage_viz_obj[accession_id]["coverage_depth"]
            }, accessions)),
            "num_accessions": taxons_to_total_accession_count[taxon]
        }

    return coverage_viz_summary_json


def get_taxons_to_total_accession_count(taxons_to_accessions):
    taxons_to_total_accession_count = {}

    for taxon, accessions in taxons_to_accessions.items():
        taxons_to_total_accession_count[taxon] = len(accessions)

    return taxons_to_total_accession_count


def select_best_accessions_per_taxon(taxons_to_accessions, accession_data, contigs_map, reads_map, num_accessions_per_taxon):
    def get_score(accession_id):
        accession_obj = accession_data[accession_id]

        contig_lengths = [contigs_map[contig_name]["alignment_length"] for contig_name in accession_obj["contigs"]]

        max_contig_length = max(contig_lengths) if len(contig_lengths) > 0 else 0
        total_contig_length = sum(contig_lengths) if len(contig_lengths) > 0 else 0
        num_reads = len(accession_obj["reads"])

        # We use the lengths of contigs, but only the number of reads, so contigs should dominate this score whenever they are present.
        return max_contig_length + total_contig_length + num_reads

    filtered_taxons_to_accessions = {}
    filtered_accession_data = {}

    for taxon, accessions in taxons_to_accessions.items():
        for accession_id in accessions:
            accession_data[accession_id]["score"] = get_score(accession_id)

        sorted_accessions = sorted(accessions, key=lambda accession_id: accession_data[accession_id]["score"], reverse=True)

        accessions_with_contigs = list(filter(lambda accession_id: len(accession_data[accession_id]["contigs"]) >= 1, sorted_accessions))

        # Take ALL accessions with one or more contigs.
        if len(accessions_with_contigs) >= num_accessions_per_taxon:
            filtered_taxons_to_accessions[taxon] = accessions_with_contigs
        # If there aren't enough accessions with contigs, take the top X.
        else:
            filtered_taxons_to_accessions[taxon] = sorted_accessions[0: num_accessions_per_taxon]

        # Move over only the accessions that were filtered. These will be processed by the pipeline step
        for accession_id in filtered_taxons_to_accessions[taxon]:
            filtered_accession_data[accession_id] = accession_data[accession_id]

    return (filtered_taxons_to_accessions, filtered_accession_data)
