#!/bin/bash

echo "\
set bgcolor white
open ranked_0.pdb
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
show #1 cartoons
view all
hide atoms
hide #1 modelsi
" >alphafold2.cxc

# Extract chain identifiers and store in an array

pdb_file="ranked_0.pdb"

# Initialize an empty array
chains=()

# Extract chain identifiers and store in an array
while IFS= read -r line; do
	chains+=("$line")
done < <(grep '^ATOM' "$pdb_file" | cut -c 22 | sort -u)


# Print the array in the desired format
echo "Chains in the PDB file: ${chains[@]}"

count=1
# Loop through the list
for (( i=0; i<${#chains[@]}; i++ )); do
	for (( j=i+1; j<${#chains[@]}; j++ )); do
		count=$((count + 1))
		echo "\
combine #1 modelId #$count name 'interface chains ${chains[i]} to ${chains[j]}'
sel #$count/${chains[i]},${chains[j]}
select ~sel & ##selected
alphafold pae #$count file ~/Desktop/pae_model_1_multimer_v3_pred_0.json plot false
alphafold contacts #$count/${chains[i]} to #$count/${chains[j]} distance 8 palette paecontacts range 0,30 radius 0.05 dashes 1
del sel
" >> alphafold2.cxc
    	done
done

