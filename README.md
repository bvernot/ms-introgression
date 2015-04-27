# ms-introgression
A short program to find introgressed haplotypes in simulated ms (or macs) data.

This requires that:
* Exactly one archaic chromosome is simulated.
* The archaic population is the last population.

## Example:

Simulate two populations - one MH and one archaic.  The two populations join at scaled time 10, with introgression between time .2 and .25.

```
ms 11 3 -s 1 -T -r 50 100000 -I 2 10 1 0 -ej 10 2 1 -em .2 1 2 3 -em .25 1 2 0 > sims.ms

python get_introgressed_regions.py -f sims.ms -jt -total -regions
```

## Options

### `-jt`
Print the time of coalescence with the archaic chromosome.  This can be used to make sure the introgressed regions can be well-separated from non-introgressed regions.

```
python get_introgressed_regions.py -f sims.ms -jt -total -regions | grep JOIN > joins.txt
```

Plot with R code:
```
library(data.table)
dt = fread('joins.txt')
hist(dt[, V2], breaks=100, col='black')
```

### -total
Print the total amount of introgressed sequence per chromosome (i.e., haplotype).

```
python get_introgressed_regions.py -f sims.ms -jt -total -regions | grep TOTAL

TOTAL pop_1 chrom_1 35977
TOTAL pop_1 chrom_2 96443
TOTAL pop_1 chrom_3 85053
TOTAL pop_1 chrom_4 65391
TOTAL pop_1 chrom_5 76507
TOTAL pop_1 chrom_6 26961
TOTAL pop_1 chrom_7 44859
TOTAL pop_1 chrom_8 26166
TOTAL pop_1 chrom_9 34165
TOTAL pop_1 chrom_10 75298
```

### -regions
Print each introgressed haplotype in bed format.  The columns are:

* **sim_pop_chrom**:
* **start**:
* **stop**:
* **chrom**:
* **pop**:
* **sim_tag**:
* **iteration_tag**:
* **intr_regions_tag**:
