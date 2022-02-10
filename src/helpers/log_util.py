
import subprocess



def log_command(logger, command, **kwargs):
    """ Return true if running the command is failed
        Return false if the command is successfully ran
    """
    logger.info(" ".join(command))

    proc = subprocess.run(" ".join(command), **kwargs)

    if proc.returncode != 0:
        logger.error("Command failed")
        if proc.stdout and isinstance(proc.stdout, str):
            logger.error(proc.stdout)
        elif isinstance(proc.stdout, bytes):
            logger.error(proc.stdout.decode())

        return True
    else:
        return False
