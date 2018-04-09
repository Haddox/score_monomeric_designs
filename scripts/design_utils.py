"""
This script contains various modules for protein design

Hugh Haddox, January 29, 2018
"""


def WriteSbatchFile(sbatch_file_name, command_file_name=None, command=None, queue_type='short', memory='2g'):
    """
    Write a file to submit a job via `sbatch`
    
    Args:
        `sbatch_file_name` : The name of the sbatch file. This will also
            serve as a prefix for the name of the output and error files.
        
        `command_file_name` : The name of a file with a list of commands
            to execute (string). The default is `None`, but if one is given,
            then this function will write an `sbatch` file that is specific
            for carying out an array of commands. If this is not given,
            the `command` argument must be specified. Only one can be
            specified at once.
            
        `command` : A command-line argument to be carried out (string).
            The default is `None`. If this is not given, the
            `command_file_name` argument must be specified. Only one can
            be specified at once.

        `queue_type` : The queue type ('short', 'medium', 'long')

	`memory` : The amount of memory in megabytes, unless other unit is specified. Default is '2g', i.e., 2 GB.
    
    Retruns:
        A file for submitting a job via `sbatch`
    """
    
    # If a file with an array of commands is provided, write an `sbatch`
    # file that is suited for this task
    if command_file_name:
        # Make sure the `command` argument hasn't been called
        assert command==None 
        
        # Write the file
        with open(sbatch_file_name, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {0}\n'.format(queue_type))
            f.write('#SBATCH --mem={0}\n'.format(memory))
            f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
            f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
            f.write('CMD=$(head -n $SLURM_ARRAY_TASK_ID {0} | tail -1)\n'.format(command_file_name))
            f.write('exec ${CMD}')
            
    # Write an `sbatch` file to carry out the specified `command`
    else:
        with open(sbatch_file_name, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {0}\n'.format(queue_type))
            f.write('#SBATCH --mem={0}\n'.format(memory))
            f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
            f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
            f.write(command)
