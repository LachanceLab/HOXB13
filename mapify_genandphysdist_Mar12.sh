awk 'BEGIN {OFS="\t"; prev_gen=0; prev_phys=0; next_gen=0; next_phys=0}
    NR==FNR {pos[$2]=$4; next}
    {
        if ($4 in pos) {
            if (prev_gen == 0) {
                prev_gen = $4; prev_phys = pos[$4];
            } else {
                next_gen = $4; next_phys = pos[$4];
                for (i = prev_gen + 1; i < next_gen; i++) {
                    if (i in interpolated) {
                        interp_phys = prev_phys + (next_phys - prev_phys) * (i - prev_gen) / (next_gen - prev_gen);
                        print interpolated[i], interp_phys, i;
                    }
                }
                delete interpolated;
                prev_gen = next_gen; prev_phys = next_phys;
            }
            print $1, pos[$4], $4;
        } else {
            interpolated[$4] = $1;
        }
    }
    END {
        for (i in interpolated) {
            print interpolated[i], prev_phys, i;
        }
    }' Path_to_pyrrohfile/YRI_recombination_map_hapmap_format_hg38_chr_22.txt pathtobimfile/chr22.bim | 
sort -k3,3n | awk '{print $1,NR,$2,$3}' > outputpath/chr22_densermapfile_genandphypos.map
