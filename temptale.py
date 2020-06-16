# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#    F I L E S   P R O C E S S I N G          #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #


import sys
import os

def check_delta_file_exist(file):
    if (os.path.exists(file)):                                                                      
        print("File " + str(file) + " is constructed.\n")
    else:
        print("There is no file " + str(file) + ". The CT-HYB needs Delta_tau.dat to run DMFT.")
        sys.exit()
        
def check_lambda_file_exist(lambda_file_name):
    if (os.path.exists(lambda_file_name)):
        print("Calculation K_tau.dat from " + str(lambda_file_name) + " with removed noise  ")
    else:
        print("There is no file ", str(lambda_file_name))
        sys.exit()

def check_ktau_file_exists(file):
    if (os.path.exists(file)):                                                                  #
        print("File " + str(file) + " is constructed.\n")                                       #
    else:                                                                                       #
        print("There is no file " + str(file) +                                                 #
        ". The CT-HYB needs Delta_tau.dat and K_tau.dat to run EDMFT.")                         #
        sys.exit()

def create_dir_with_files(type_of_calc, iteration):
    folder = "results_iteration_" + str(iteration)
    os.system("mkdir " + folder)
    os.system("cp ./*.dat " + folder )
    os.system("cp ./*.pdf " + folder )
    os.system("cp ./*.h5 " + folder )
    os.system("cp ./params " + folder )
    os.system("rm -fr *.dat")
    os.system("rm -fr *.pdf")
    os.system("rm -fr *.h5")
    os.system("rm params")

def prepare_files_for_new_it(type_of_calc, iteration):
    folder = "results_iteration_" + str(iteration)
    os.system("cp " + folder + "/Delta_new_minimized.dat Delta.dat")
    
    if(type_of_calc == "edmft"):
        os.system("cp " + folder + "/Lambda_new_smooth.dat Phi.dat")
