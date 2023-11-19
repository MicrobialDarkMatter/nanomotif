import logging as log

def configure_logger(log_file, verbose=False):
    if verbose:
        level=log.DEBUG
    else:
        level=log.INFO
    
    log.basicConfig(
        filename=log_file,
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )

