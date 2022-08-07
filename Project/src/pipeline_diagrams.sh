# This produces diagrams of all the steps and files in the project's pipeline,
# regardless of whether those steps need to be run or whether the files exist yet.
# Run this from the project directory ('sh src/diagrams.sh')
snakemake --rulegraph | dot -Tpng > rulegraph.png
snakemake --filegraph | dot -Tpng > filegraph.png
