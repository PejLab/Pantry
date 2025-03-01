# This produces diagrams of all the steps and files in the project's pipeline,
# regardless of whether those steps need to be run or whether the files exist yet.
# Run this from the project directory ('sh scripts/pipeline_diagrams.sh')
snakemake --rulegraph | dot -Tpng > diagram_steps.png
# snakemake --filegraph | dot -Tpng > diagram_files.png
