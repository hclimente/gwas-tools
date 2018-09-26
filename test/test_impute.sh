impute --ped data/example.ped --map data/example.map --reference data/1000GP_Phase3 -resume --strand_info data/strand_info --population EUR "$@"
rm -f out.ped out.map