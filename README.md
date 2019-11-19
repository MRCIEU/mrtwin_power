# Notes from call on 5/12/2018


## Traits to analyse:

```
height -> education
bmi -> education
height -> high blood pressure
bmi -> high blood pressure
height -> type 2 diabetes
bmi -> type 2 diabetes
```

In HUNT Ben has done `ldl -> chd` but this is not available in UKBB

Ben mentioned that there might be a set of examples that are more indicative of dynastic effects. e.g. 

```
Alcohol -> education?
Parent's age of death -> body mass index?
```


## Methods to run

For the set of traits we need to decide on how to implement the methods.



### OPTION A - more robust to horizontal pleiotropy

Every SNP is gives an estimate using:

- within family - PLM method with robust SE
- population (one sample per family) - Standard IV

Then we meta analyse to get

- Overall PLM
- Overall IV



### OPTION B - more robust to weak instrument bias and better powered

Generate risk score for the exposure as instrument and run:

- within family - PLM method with robust SE
- population (one sample per family) - Standard IV

We can then discuss that if weak instrument bias was not an issue we could extend to the modular approach of OPTION A

If there is weak instrument bias here it will be towards observational estimate



### OPTION C - more robust to horizontal pleiotropy and sample overlap issues

Do the following

1. Get SNP-exposure effects on height in HUNT for 385 instruments
2. Get SNP-outcome effects on education in UKBB for instruments in (1)
3. Harmonise
4. MR IVW + Mode + Median etc


1 and 2 can be done
- Using just one individual per family
- Using basic PLM (no instrument) or DiD
	- `plm(height ~ rs123, index="familyID", model="within")`

Ben (HUNT) and Neil (UKBB) to generate file of every SNP against every trait including

```
- trait
- SNP
- beta
- se
- pvalue 
- effect allele
- non-effect allele
- effect allele frequencies
- sample size
```

If there is weak instrument bias here it will be towards the null


### OPTION D - Alternative to C for the single sample context

Single sample setting, assumes strong instruments (?)

For each SNP:

```
plm(outcome ~ exposure | inst, index="famid")
```

This will give the causal effect estimate based on each SNP in a within family setting. If there are p instruments then we can combine those p estimates using small adaptations to standard two-sample MR pleiotropy tools e.g.


- Radial MR approach:
	- IVW fixed  I(b_wr * sqrt(se)) ~ 0 + sqrt(se)
	- IVW random
	- heterogeneity
	- Egger
- Median
- Mode

Weak instruments will bias towards observational association but the method will account for sample overlap problems



### OPTION E - improve on power in D

Same as D but do it in both HUNT and UKBB, then meta analyse each SNP first so

`x -> y` estimate using rs123 will be made in UKBB and HUNT, so then meta analyse to get an overall estimate. Then use those in the framework as in option D



### To run

```
module add languages/anaconda3/5.2.0-tflow-1.11
snakemake -prk \
-j 400 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```

