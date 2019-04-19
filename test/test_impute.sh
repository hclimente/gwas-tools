impute --bfile data/example --reference ${HOME}/data/1000GP_Phase3 -resume --strand_info data/strand_info.txt --population EUR "$@"
rm -f out.ped out.map
