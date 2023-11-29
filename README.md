# alphafold

This is a script to submit AlphaFold2-Multimer jobs with minimal user input in a SLURM compute environment.
Once installed, add executable permissions with

chmod +x alphafold/alphafold.sh

Then,

./alphafold/alphafold.sh [options]

Options:
  -d DIRECTORY   Specify the directory path. (Required)
  -f FILE        Specify the sequence file path. (Required)
  -m MODE        Specify 'quick' or 'thorough'. ([1] or 5 predictions/model).
  -h             Show help information.

You could add a soft link to the script in your home folder, for example, by

ln -s alphafold/alphafold.sh alphafold

The script will assess if the sequence file provided is in FASTA format and contains only protein sequences.
The sequence file is provided with -f flag, for example, -f fastafile.fa

The script will also determine if the directory is writeable, and the output will be in a subdirectory with the user's login name.

The -m flag for a quick or thorough mode of prediction: 
- 'quick' generates 1 prediction per model (5 total)
- 'thorough' generates 5 predictions per model (25 total)

There will be a basic determination of required compute resources and time limits based on the length of the protein sequences:

<249 residues       6 hr 
250-999 residues    24 hr 
1000-1299 residues  36 hr
1300-2499 residues  48 hr
>2500 residues      72 hr

For predictions >2500 residues, only the 'quick' mode of 1 prediction/model is permitted, as these are typically very long jobs.

Once structure predictions are available, the script will run ChimeraX and generate a short movie of the top prediction.
ImageMagick's 'convert' command will convert this to a GIF. 

A ChimeraX script with the extension _chimera_align.cxc will align all predictions on the top solution with the largest chain as the reference.
A predicted aligned error (PAE) plot will also be generated.

Another ChimeraX script with the extension _chimera_pae.cxc will generate pseudobonds between residues of each pair of chains in the top prediction.
The pseudobonds will be colored on a spectrum of PAE values (blue-red is 0-30).

Finally, the script will zip all the predictions, the JSON file with PAE scores for the top one, and the two ChimeraX scripts above and email them to the user.
The animated GIF with the rotating top prediction will also be inserted into the body of the email.

Afterward, if the script picks up Box.com credentials in the .netrc, the script will directly upload the entire structure prediction to Box for long-term storage.
Note that these directories can be significant in size, especially for large predictions.

Any questions or problems with the script can be reported to the author. Email can be found by searching for the author's name and affiliation below.

Eugene Valkov, D.Phil.
Center for Cancer Research
National Cancer Institute
National Institutes of Health
Frederick, Maryland, U.S.A.
