import logging
import shlex
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def command_text(command):
    return shlex.join(str(part) for part in command)


def run_command(command, dry_run=False, stdout=None, stderr=None):
    command = [str(part) for part in command]
    logging.info("Command: %s", command_text(command))
    if dry_run:
        return None
    return subprocess.run(command, check=True, stdout=stdout, stderr=stderr, text=False)


def require_executable(name):
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required executable not found in PATH: {name}")
    return path


def require_file(path, label="file"):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Required {label} not found: {path}")
    return path


def require_bowtie_index(prefix):
    prefix = Path(prefix)
    suffixes = (".1.ebwt", ".1.ebwtl")
    if not any(Path(str(prefix) + suffix).is_file() for suffix in suffixes):
        raise FileNotFoundError(f"Bowtie index prefix not found: {prefix}")


def sample_name(path, suffix=None):
    name = Path(path).name
    candidates = (suffix,) if suffix else FASTQ_SUFFIXES
    for candidate in candidates:
        if candidate and name.endswith(candidate):
            name = name[:-len(candidate)]
            break
    return name.replace("-", "_")


def discover_fastqs(input_dir, suffix):
    input_dir = Path(input_dir)
    files = sorted(path for path in input_dir.glob(f"*{suffix}") if path.is_file())
    if not files:
        raise FileNotFoundError(f"No '*{suffix}' files found in {input_dir}")
    samples = [(sample_name(path, suffix), path) for path in files]
    names = [name for name, _ in samples]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate normalized sample names: {', '.join(duplicates)}")
    return samples


def run_parallel(worker, items, config, jobs, label):
    if not items:
        logging.info("No samples selected for %s", label)
        return
    if config.get("dry_run") or jobs == 1:
        for item in items:
            worker(item, config)
        return

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(worker, item, config): item[0] for item in items}
        for future in as_completed(futures):
            sample = futures[future]
            try:
                future.result()
            except Exception as error:
                for pending in futures:
                    pending.cancel()
                raise RuntimeError(f"{label} failed for sample '{sample}'") from error
