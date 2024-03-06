mkdir output
mkdir output/dada2_pipeline
mv quality_forwards.pdf output
mv quality_reverse.pdf output
mv error_model_forwards.pdf output
mv error_model_reverse.pdf output
mv seqtab.rds output
mv seqtab_nochim.rds output
mv sequence_process_summary.txt output
mv taxID.rds output
cp -r ../../dada2_pipeline/ output/dada2_pipeline