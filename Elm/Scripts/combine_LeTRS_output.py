"""
Combine all files of a certain type into a single output file. 
This is a custom script to LeTRS output files, adding a column for the source file.
The overall logic should work for other file types as well with minor modifications.

Written by: Elm M.
"""
import os
import regex as re


def find_file_types(pattern):
    matched_files = []
    for root, dirs, files in os.walk('.'):
        for file in files:
            if re.search(pattern, file):
                matched_files.append(os.path.join(root, file))
    return matched_files

def combine_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for f in range(len(file_list)):
            with open(file_list[f], 'r') as infile:
                infile_lines = infile.readlines()
                infile.close()
                for l in range(len(infile_lines)):
                    if l == 0:
                        if f == 0:
                            # Add column for which source file the line is from
                            infile_lines[l] = infile_lines[l].rstrip() + "\tSource_File\n"
                            outfile.write(infile_lines[l])
                            continue
                        continue
                    # Specific to LeTRS
                    elif infile_lines[l].startswith('The'):
                        break
                    else:
                        file_name = re.search(r'(?<=\.\/).+(?=\/results)', file_list[f]).group(0).replace(r'/', '_').replace('.', '_')
                        infile_lines[l] = infile_lines[l].rstrip() + f"\t{file_name}\n"
                        outfile.write(infile_lines[l])
    outfile.close()


def main():
    # Set working directory
    os.chdir('/restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/LeTRS_Output/')
    # File types you want to combine
    patterns = {"cannonical": 'known_junction.tab',\
                "nc-Leader": 'novel_junction.tab',\
                "nc-ind": 'TRS_L_independent_junction.tab'}
    for key, pattern in patterns.items():
        files_to_combine = find_file_types(pattern)
        output_filename = f'/restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/LeTRS_Output/combined/combined_{key}.tab'
        combine_files(files_to_combine, output_filename)

if __name__ == "__main__":
    main()