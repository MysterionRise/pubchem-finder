import argparse
import gzip
import logging
import os
import pathlib
from datetime import datetime
from typing import Callable, Dict, Generator, List

import ftpretty

from core import FTP, info

STATE_DIR = 'state'


class Extractor:
    def __init__(self, chembl_dir):
        self.chembl_dir = pathlib.Path(chembl_dir)

    @staticmethod
    def __gunzip_content(file: pathlib.Path) -> pathlib.Path:
        target_file = file.with_suffix('')
        with open(target_file, 'wb') as target, gzip.open(file, 'rb') as src:
            target.write(src.read())
        return target_file

    @staticmethod
    def cast_path(
        root: pathlib.Path, files: List[str]
    ) -> Generator[pathlib.Path, None, None]:
        for file in files:
            yield root / file

    def extract(self, handler: Callable[[pathlib.Path], None]):
        for root, dirs, files in os.walk(self.chembl_dir):
            for file in self.cast_path(pathlib.Path(root), files):
                if file.suffix == '.gz':
                    try:
                        info(f'Extracting file {file}')
                        extracted = self.__gunzip_content(file)
                        handler(extracted)
                        extracted.unlink(missing_ok=True)
                    except EOFError as err_:
                        logging.error('File %s: %s', file, err_)
                    finally:
                        info(f'File {file} uploaded')


class State:
    def __init__(self, file: pathlib.Path):
        self.file = file
        try:
            with open(self.file, 'x'):
                pass
        except FileExistsError:
            pass

    def get_state(self):
        pass

    def serialize_state(self, state: str) -> Dict[str, str]:
        if not state:
            return {'state': '', 'last_change': '', 'message': ''}
        exploded = state.split('\n')
        return {
            'state': exploded[0],
            'last_change': exploded[1],
            'message': exploded[2],
        }

    def unserialize_state(self, state: Dict[str, str]) -> str:
        result = ''
        for value in state.values():
            result += f'{value}\n'
        return result

    def change_state(self, state: str, message: str = ''):
        with open(self.file, 'r+') as state_file:
            cur_state = self.serialize_state(state_file.read())
            cur_state['state'] = state
            cur_state['message'] = message
            cur_state['last_change'] = datetime.utcnow()
            state_file.seek(0)
            state_file.truncate()
            str_state = self.unserialize_state(cur_state)
            state_file.write(str_state)


def _extract(
    filename: pathlib.Path,
    tmpdir: pathlib.Path,
    handler: Callable[[pathlib.Path], None],
):
    source_file = tmpdir / filename
    if source_file.suffix == '.gz':
        target_file = source_file.with_suffix('')
        try:
            with open(target_file, 'wb') as target, gzip.open(
                source_file, 'rb'
            ) as src:
                target.write(src.read())
            handler(target_file)
        finally:
            target_file.unlink(missing_ok=True)


def _download(
    ftp: ftpretty,
    filename: pathlib.Path,
    tmpdir: pathlib.Path,
    max_tries: int = 20,
    tries: int = 0,
):
    if tries == max_tries - 1:
        raise TimeoutError('Max tries exceeded')


def _checksum(target_file: pathlib.Path, tmpdir: pathlib.Path) -> bool:
    return True


class Pull:
    def __init__(
        self,
        state: Dict,
        workdir: pathlib.Path,
        tmpdir: pathlib.Path,
        timeout: int,
    ):
        self.state = state
        self.workdir = workdir
        self.tmpdir = tmpdir
        self.timeout = timeout
        self.ftp = 'ftp.ncbi.nlm.nih.gov'

    def execute(self):
        source_dir = (
            pathlib.Path('pubchem') / 'Compound' / 'CURRENT-Full' / 'SDF'
        )
        with FTP(
            self.ftp, 'anonymous', 'anonymous@domain.com', timeout=self.timeout
        ) as ftp_:
            for file in sorted(ftp_.list(str(source_dir))):

                file = pathlib.Path(file)

                if file.suffix != '.gz':
                    continue

                info(f'Downloading {file.name}')

                state_file = (
                    self.workdir
                    / STATE_DIR
                    / pathlib.Path(file.name).with_suffix('.state')
                )

                try:
                    _download(ftp_, file, self.tmpdir)
                except TimeoutError:
                    State(state_file).change_state('timeouterror')
                else:
                    State(state_file).change_state('download_ok')

                info('Ok')

                info(f'Checking {file.name}')

                if not _checksum(file, self.tmpdir):
                    info('ERROR')
                    State(state_file).change_state('checksum_error')
                    continue
                State(state_file).change_state('checksum_ok')
                info('Ok')

                info(f'Extracting {file.name} into Elastic')

                _extract(
                    pathlib.Path(file.name),
                    self.tmpdir,
                    lambda x: print(f'Hello {x}'),
                )
                State(state_file).change_state('extract_ok')
                info('Ok')


def create_dirs(workdir: pathlib.Path, tmpdir: pathlib.Path):
    dirs_ = [
        workdir,
        workdir / STATE_DIR,
        tmpdir,
    ]
    for dir_ in dirs_:
        dir_.mkdir(parents=True, exist_ok=True)


def read_state(workdir: pathlib.Path) -> Dict[str, str]:
    statedir = workdir / STATE_DIR
    files = [pathlib.Path(name).name for name in statedir.glob('*.state')]
    state = {}
    for file in files:
        state[file] = file

    return state


def pull(args: argparse.Namespace) -> None:
    workdir = pathlib.Path(args.workdir)
    tmpdir = pathlib.Path(args.tmpdir)
    create_dirs(workdir, tmpdir)
    puller = Pull(read_state(workdir), workdir, tmpdir, 30)
    puller.execute()
