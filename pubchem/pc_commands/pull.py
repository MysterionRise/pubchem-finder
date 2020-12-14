import argparse
import ftplib
import gzip
import json
import logging
import pathlib
import time
from datetime import datetime
from typing import Callable, Dict, Generator, List

from core import FTP, info
from elastic import ElasticDatabase

STATE_DIR = "state"


class State:
    # TODO: move states from files to elastic
    def __init__(self, file: pathlib.Path, source_file: pathlib.Path = None):
        self.file = file
        try:
            with open(self.file, "x") as state_file:
                if not source_file:
                    raise ValueError(
                        "source file cannot be none, "
                        "when state file is not exists"
                    )
                self.__write_state(
                    state_file,
                    {"source_file": str(source_file), "state": "download"},
                )
        except FileExistsError:
            pass

    def get_state(self) -> Dict:
        with open(self.file, "r") as state_file:
            return self.__load(state_file.read())

    def __load(self, state: str) -> Dict[str, str]:
        return json.loads(state)

    def change_state(self, state: str, message: str = ""):
        with open(self.file, "r+") as state_file:
            cur_state = json.load(state_file)
            cur_state["state"] = state
            cur_state["message"] = message
            cur_state["last_change"] = datetime.now().strftime(
                "%Y-%m-%d %H:%M:%S"
            )
            state_file.seek(0)
            state_file.truncate()
            self.__write_state(state_file, cur_state)

    def __write_state(self, state_file, state: Dict):
        json.dump(state, state_file)


def _extract(
    filename: pathlib.Path,
    tmpdir: pathlib.Path,
    handler: Callable[[pathlib.Path], None],
):
    source_file = tmpdir / filename
    if source_file.suffix == ".gz":
        target_file = source_file.with_suffix("")
        try:
            with open(target_file, "wb") as target, gzip.open(
                source_file, "rb"
            ) as src:
                target.write(src.read())
            handler(target_file)
        finally:
            target_file.unlink(missing_ok=True)


class Pull:
    def __init__(
        self,
        state: Dict,
        workdir: pathlib.Path,
        tmpdir: pathlib.Path,
        elastic: ElasticDatabase,
        timeout: int,
    ):
        self.state = state
        self.workdir = workdir
        self.tmpdir = tmpdir
        self.timeout = timeout
        self.ftp = "ftp.ncbi.nlm.nih.gov"
        self.elastic = elastic

        self.ftp_files = []

    def __get_states(self) -> Generator[State, None, None]:
        file: pathlib.Path
        for file in sorted(self.ftp_files):
            state_file = self.state.get(file.name)
            if not state_file:
                state_file = (
                    self.workdir / STATE_DIR / file.with_suffix(".state").name
                )
            yield State(state_file, file)

    def __download(self, state: State, tries=0) -> str:
        source_file = pathlib.Path(state.get_state().get("source_file"))
        local_file = self.tmpdir / source_file.name
        if tries >= 20:
            state.change_state("download_error", "Max tries exceeded")
            raise TimeoutError("Max tries exceeded")
        with FTP(
            self.ftp, "anonymous", "anonymous@domain.com", timeout=self.timeout
        ) as ftp_:
            try:
                with open(local_file, "xb") as target_file:
                    info(f"Download {source_file}")
                    ftp_.get(source_file, target_file)
            except FileExistsError:
                local_file.unlink()
                self.__download(state, tries + 1)
            except (ftplib.error_temp, EOFError, BaseException) as err_:
                logging.warning("%s, sleeping 5 seconds", err_)
                time.sleep(5)
                logging.warning("retry...")
                self.__download(state, tries + 1)
        next_state = "checksum"
        state.change_state(next_state, "")
        info(f"Ok")
        return next_state

    def __checksum(self, state: State):
        source_file = pathlib.Path(
            state.get_state().get("source_file") + ".md5"
        )
        local_file = self.tmpdir / source_file.name

        info(f"Check md5 {source_file}")
        with FTP(
            self.ftp, "anonymous", "anonymous@domain.com", timeout=self.timeout
        ) as ftp_:
            with open(local_file, "wb") as target_file:
                info(f"Download {source_file}")
                ftp_.get(source_file, target_file)
        reference_sum = ""
        with open(local_file, "r") as md5_ref_file:
            reference_sum = md5_ref_file.read().split(" ")[0]

        # TODO: check md5 in subprocess

        logging.warning("Skip md5 checksum")

        next_state = "extract"
        state.change_state(next_state, "")
        return next_state

    def __extract(self, state: State):
        source_file = (
            self.tmpdir
            / pathlib.Path(state.get_state().get("source_file")).name
        )
        target_file = self.tmpdir / source_file.with_suffix("").name
        info(f"Extract {source_file}")

        with open(target_file, "wb") as target, gzip.open(
            source_file, "rb"
        ) as src:
            target.write(src.read())

        next_state = "load"
        state.change_state(next_state, "")
        info(f"Ok")
        return next_state

    def __load(self, state: State):
        gzipped_file = (
            self.tmpdir
            / pathlib.Path(state.get_state().get("source_file")).name
        )
        text_file = self.tmpdir / gzipped_file.with_suffix("").name
        info(f"Load {text_file}")

        self.elastic.handler(text_file)

        next_state = "clear"
        state.change_state(next_state, "")
        info(f"Ok")
        return next_state

    def __clear(self, state: State):
        gzipped_file = (
            self.tmpdir
            / pathlib.Path(state.get_state().get("source_file")).name
        )
        text_file = self.tmpdir / gzipped_file.with_suffix("").name
        target_md5_file = pathlib.Path(
            state.get_state().get("source_file") + ".md5"
        )
        local_md5_file = self.tmpdir / target_md5_file.name
        info(f"Clear")
        gzipped_file.unlink(missing_ok=True)
        text_file.unlink(missing_ok=True)
        local_md5_file.unlink(missing_ok=True)
        next_state = "complete"
        state.change_state(next_state, "")
        info(f"Ok")
        return next_state

    def __run_from_state(self, state: State):

        cur_state = state.get_state().get("state")
        if cur_state in ("download", "download_error", "checksum_error"):
            cur_state = self.__download(state)
        if cur_state == "checksum":
            cur_state = self.__checksum(state)
        if cur_state == "extract":
            cur_state = self.__extract(state)
        if cur_state == "load":
            cur_state = self.__load(state)
        if cur_state == "clear":
            # todo: clear
            self.__clear(state)

    def execute(self):
        source_dir = (
            pathlib.Path("pubchem") / "Compound" / "CURRENT-Full" / "SDF"
        )

        # Get all file names
        with FTP(
            self.ftp, "anonymous", "anonymous@domain.com", timeout=self.timeout
        ) as ftp_:
            for file in ftp_.list(str(source_dir)):
                file = pathlib.Path(file)
                if file.suffix == ".gz":
                    self.ftp_files.append(file)

        # Start job for every state
        for state in self.__get_states():
            # todo: add multithreading
            try:
                self.__run_from_state(state)
            except Exception as err_:
                logging.error(
                    "Error %s. Additional info in %s", err_, state.file
                )


def create_dirs(workdir: pathlib.Path, tmpdir: pathlib.Path):
    dirs_ = [
        workdir,
        workdir / STATE_DIR,
        tmpdir,
    ]
    for dir_ in dirs_:
        dir_.mkdir(parents=True, exist_ok=True)


def read_state(workdir: pathlib.Path) -> Dict[str, pathlib.Path]:
    """Generates pairs file.gz : path / to / state_file.state
    :param workdir:
    :return:
    """
    state_dir = workdir / STATE_DIR
    files = [pathlib.Path(name).name for name in state_dir.glob("*.state")]
    state = {}
    for file in files:
        state[pathlib.Path(file).with_suffix(".gz").name] = state_dir / file

    return state


def pull(args: argparse.Namespace) -> None:
    workdir = pathlib.Path(args.workdir)
    tmpdir = pathlib.Path(args.tmpdir)
    create_dirs(workdir, tmpdir)

    puller = Pull(
        read_state(workdir), workdir, tmpdir, ElasticDatabase(args), 30
    )
    puller.execute()
