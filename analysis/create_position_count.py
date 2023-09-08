import gzip
import json
import sys

def get_header_positions(filename: str) -> tuple[int, int]:

    with gzip.open(filename, mode="rt") as file:
        header = file.readline()

    # get position of normal_genotype and tumor_genotype
    header_list = header.split("\t")

    try:
        normal_gt_index = header_list.index("normal_genotype")
        tumor_gt_index = header_list.index("tumor_genotype")
    except ValueError:
        print("element not in header")

    return (normal_gt_index, tumor_gt_index)


def count_occurences(filename: str, normal_gt_index: int, tumor_gt_index) -> dict:

    position_count = {}

    with gzip.open(filename, mode="rt") as file:
        for i, line in enumerate(file):
            #if i > 100000000:
            #    break

            sys.stdout.write("Progress: working on line %d   \r" % (i) )
            sys.stdout.flush()

            splitline = line.split("\t")

            normal_gt = splitline[normal_gt_index]
            tumor_gt = splitline[tumor_gt_index]

            if normal_gt == "1,1" and (tumor_gt.startswith("0,") or ",0" in tumor_gt) and tumor_gt != "0,0":
                key = splitline[0] + "-" + splitline[1]
                if key in position_count.keys():
                    position_count[key] += 1
                else:
                    position_count[key] = 1
    
    return position_count

def main(filename: str, outfile: str):

    normal_gt_index, tumor_gt_index = get_header_positions(filename)
    position_count_dict = count_occurences(filename, normal_gt_index, tumor_gt_index)

    with open(outfile, "w") as file:
        json.dump(position_count_dict, file)

if __name__ == "__main__":

    filename = "mmml_all_patients_extra_columns.vcf.gz"
    outfile = "LOH_position_count_mmml_all_patients_extra_columns.json"

    main(filename, outfile)
    