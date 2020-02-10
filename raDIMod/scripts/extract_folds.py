file_fold_name = ["domains_found.txt","families_found.txt","folds_found.txt","superfamilies_found.txt"]


for file_name in file_fold_name:
    file = open(file_name).readlines()
    for lines in file:
        k = lines.split()
        print(k)
