#renaming mouse transcriptomics files
import os

path_folder = "/restricted/projectnb/challenge2025/Data/mouse_transcriptomics"

sample_num = 0

for file_name in os.listdir(path_folder):
    """
    current_path = os.path.join(path_folder, file_name)
    split_name = file_name.split("_")
    num = (split_name[0])[-2:]
    split_ind = (split_name[1]).split(".")[0]
    #print(num)

    new_name = split_name[0] + "-r" + split_ind + ".fq.gz"
    new_path = os.path.join(path_folder, new_name)
    print(new_path)
    #os.rename(current_path, new_path)
    """
    sample_num+=1
print(sample_num)