def get_server_run():
    path_to_exec_file   = '/storage/praha1/home/medvedeva/workspace/other_codes/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads     = 16
    path_to_maxent      = ''
    return path_to_exec_file, num_mpi_threads, path_to_maxent

def get_local_run():
    path_to_exec_file   = '/Users/witcher/workspace/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads     = 3
    path_to_maxent      = '/Users/witcher/workspace/CT_HYB_SEGMENT/Maxent/build2/maxent'
    return path_to_exec_file, num_mpi_threads, path_to_maxent
