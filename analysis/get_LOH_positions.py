import gzip
import sys

def write_LOH_positions(all_patients_file: str, outfile: str):

    with gzip.open(all_patients_file, mode="rt") as read_file:
        with gzip.open(outfile, mode="wt") as print_file:
            for i, line in enumerate(read_file):
                
                sys.stdout.write("Progress: working on line %d   \r" % (i) )
                sys.stdout.flush()

                splitline = line.split("\t")

                if i == 0:
                    header_list = splitline

                    try:
                        normal_gt_index = header_list.index("normal_genotype")
                        tumor_gt_index = header_list.index("tumor_genotype")

                    except ValueError:
                        print("element not in header")
                        break
                
                normal_gt = splitline[normal_gt_index]
                tumor_gt = splitline[tumor_gt_index]

                if normal_gt == "1,1" and (tumor_gt.startswith("0,") or ",0" in tumor_gt) and tumor_gt != "0,0":
                    print_file.write(line)

def main(all_patients_file: str, outfile: str):
    write_LOH_positions(all_patients_file, outfile)

if __name__ == "__main__":
    all_patients_file = "mmml_all_patients_extra_columns.vcf.gz"
    outfile = "LOH_positions_mmml.vcf.gz"
    main(all_patients_file, outfile)