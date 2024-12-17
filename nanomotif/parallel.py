import multiprocessing
import progressbar
import time


def update_progress_bar(progress, total_tasks, print_progress, timeout=30000):
    """
    Update progress bar for parallel processing

    Parameters:
    - progress (multiprocessing.Value): The progress counter
    - total_tasks (int): The total number of tasks
    - timeout (int): The timeout in seconds
    """
    last_update_time = time.time()
    with progressbar.ProgressBar(max_value=total_tasks) as bar:
        while True:
            current_time = time.time()
            if progress.value >= total_tasks:
                break
            if current_time - last_update_time > timeout:
                print("Timeout reached. Exiting progress bar update.")
                break
            if print_progress:
                bar.update(progress.value)
            time.sleep(0.5)
