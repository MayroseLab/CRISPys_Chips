import globals

def createHeaderJob(path, job_name, queue="itaym", ncpu=1, mem=16):
    """
    A function to create qsub file with activating crispys conda env and ssh to 0-247 machine

    :param path: path to log files
    :param job_name: job name
    :param queue:
    :param ncpu: cpu number default 1
    :param mem: memory to use (in gb) default 16
    :return: a string that can be used to write sh file to run on the cluster (need to add command before running on the cluster)
    """
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += f"#PBS -q {queue}\n"
    text += "#PBS -N " + job_name + "\n"
    text += "#PBS -e " + path + "/" + job_name + ".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name + ".OU" + "\n"
    text += "#PBS -l select=ncpus=" + str(ncpu) + ":mem=" + str(mem) + "gb\n"
    text += "source ~/.bashrc\n"
    text += "export PATH='$CONDA_PREFIX/bin:$PATH'\n"
    text += f"conda activate {globals.conda_environment}"
    return text