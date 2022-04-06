import pandas as pd
import pymol
from pymol import cmd
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("inputDir", type=str, help="directory with .pdb or .cif files")
parser.add_argument("outputDir", type=str, help="name of the directory to save results to (structures and dataset)")
parser.add_argument("threshold", type=int, default="15", help="maximum length (in residues) of biological polymer to keep, default = 15")
parser.add_argument("-p", "--peptide", action="store_true", help="if true, only keep biological polymers - delete small molecules, solvents, and ions")
parser.add_argument("-c", "--clean", action="store_true", help="if true, clean up files by removing extra chains")
parser.add_argument("-b", "--bsa", action="store_true", help="if true, calculate SASA and BSA")
args = parser.parse_args()

cols = ['pdb_id', 'chain 1', 'chain 2', 'protein_sasa', 'peptide_sasa', 'combined_sasa', 'bsa']
lst = []
path = args.inputDir
dirName = args.outputDir

# Create target Directory if don't exist
if not os.path.exists(dirName):
    os.mkdir(dirName)
else:    
    print("Directory " , dirName ,  " already exists")


#######################################################################################################################################

def pdb_separator():
    
    # loop through files and calculates SASAs and BSA #
    for filename in os.listdir(path):
        if filename.endswith(".pdb") or filename.endswith(".cif"):
            filenamefull = os.path.join(path, filename)
            chain_lens = []
            chain_lens_ca = []
            i = 0
            repetition_check = 0
            try:
                # reloads pymol and loads file #
                cmd.reinitialize()
                cmd.load(filenamefull)

                # gets PDB ID without extension #
                file = os.path.splitext(filename)

                # gets all chains #
                chains = cmd.get_chains(file[0])
                
                print("Cleaning chains for file ", file[0])

                # looks for the chain with most atoms (assumes this is the most complete chain of the protein) #
                # then, delete all other chains that have less than X alpha carbons (thus, keeps peptides)    #
                for chain in chains:
                    count = cmd.count_atoms(file[0] + ' and chain ' + chain)
                    count_ca = cmd.count_atoms(file[0] + ' and chain ' + chain + ' and name CA')
                    chain_lens.append(count)
                    chain_lens_ca.append(count_ca)

                for n in chain_lens:
                    if n < max(chain_lens) and chain_lens_ca[i] > args.threshold:
                        cmd.remove("chain "+chains[i])
                        i = i + 1
                    elif n == max(chain_lens) and repetition_check == 0:
                        repetition_check = 1
                        i = i + 1
                    elif n == max(chain_lens) and repetition_check == 1:
                        cmd.remove("chain "+chains[i])
                        i = i + 1
                    else:
                        i = i + 1

                max_value = max(chain_lens)
                max_index = chain_lens.index(max_value)

                # gets the peptide that is interacting with the leftover chain, selects the complex and deletes the rest #
                if args.peptide:
                    cmd.create("complex", (file[0] + ' and chain '+chains[max_index]+' + bychain all within 3 of chain '+chains[max_index]))
                    cmd.create("complex", "complex and polymer")
                else:
                    cmd.create("complex", (file[0] + ' and chain '+chains[max_index]+' + bymolecule all within 3 of chain '+chains[max_index]))
                cmd.remove("not complex")
                filename = 'clean_' + filename
                cmd.save(os.path.join(dirName, filename),"all",-1)
            
            except:
                print("failed at cleaning PDB ID", filename) 
    return(0)

########################################################################################################################################################         

def calculate_bsa():
    
    if not args.clean:
	    path2 = path
    else:
	    path2 = os.path.join(path, dirName)

    for filename in os.listdir(path2):

        if filename.endswith(".pdb") or filename.endswith(".cif"):
            print(filename)
            try:
                filenamefull = os.path.join(path2, filename)
                cmd.reinitialize()
                cmd.load(filenamefull)

                # gets PDB ID without extension #
                file = os.path.splitext(filename)
                # gets first two chains after removing duplicate entries #
                chains = cmd.get_chains(file[0])

                # create objects for protein, peptide, and protein+peptide pair #
                cmd.create("protein", (file[0] + ' and chain '+chains[0]))
                cmd.create("peptide", (file[0] + ' and chain '+chains[1]))
                cmd.create("combined", (file[0] + ' and chain '+chains[0]+'+'+chains[1]))

                # add missing hydrogens #
                cmd.h_add()

                # ignores solvent for calculations #
                cmd.flag("ignore", "none")
                cmd.flag("ignore", "solvent")

                # sets SASA calculation parameters #
                cmd.set("dot_solvent", 1)
                cmd.set("dot_density", 3)
                
                print("Calculating areas for file ", file[0])

                # calculates areas #
                protein_area=cmd.get_area("protein")
                peptide_area=cmd.get_area("peptide")
                combined_area=cmd.get_area("combined")

                # add to list #
                lst.append([file[0], chains[0], chains[1], protein_area, peptide_area, combined_area, ((protein_area + peptide_area) - combined_area)])

            except:
                print("failed at PDB ID", filename) 
    
    df1 = pd.DataFrame(lst, columns=cols)
    dataset_name = "areas.csv"
    df1.to_csv(os.path.join(dirName, dataset_name))
                
########################################################################################################################################################

if args.clean:
    pdb_separator()

if args.bsa:
    calculate_bsa()
    
if not args.bsa and not args.clean:
    print("You need to pass at least one extra argument, -c to clean extra chains and/or -b to calculate SASA and BSA")