#Genome annotation by gutSMASH
python  run_gutsmash.py --cb-knownclusters --enable-genefunctions  data_path/genome_file  --output-dir  examples/caiTABCDE/results_genome_file_name  --genefinding-tool  prodigal  -c  16

#Sequences extraction
python3  Gut-extractor.py  -n  carnitine_degradaion_caiTABCDE  -i  examples/caiTABCDE  -o  caiTABCDE
