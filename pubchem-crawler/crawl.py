import ftplib
import logging
import time
from pathlib import Path
from typing import Generator, Iterable, List, Set
import argparse

from ftpretty import ftpretty


class FTP:
    def __init__(self, *args, **kwargs):
        self.conn = ftpretty(*args, **kwargs)

    def __enter__(self):
        return self.conn

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.conn.close()
        except Exception as e:
            logging.error(e)


class PubchemLoader:

    dir_mapping = {"compounds": "Compound"}

    # Order make sense
    sup_extensions = [".md5", ".gz"]

    def __init__(self, data: str, target_dir: str, fmt: str) -> None:
        self.ftp = "ftp.ncbi.nlm.nih.gov"
        self.target = Path("pubchem") / self.dir_mapping[data]
        self.fmt = fmt.upper()
        self.target_dir = Path(target_dir)

    def __checksum(self, target_file: Path) -> bool:
        # todo: create checksum functionality
        target_file.with_suffix(".gz.md5").is_file()
        return True

    def __download_file(self, ftp_: ftpretty, filename: Path, tries: int = 0):
        pth = self.target_dir / filename
        if pth.suffix not in self.sup_extensions:
            return
        try:
            pth.parent.mkdir(parents=True)
        except FileExistsError:
            pass
        try:
            with open(pth, "xb") as target_file:
                ftp_.get(filename, target_file)
                logging.info("Downloaded %s", pth)
        except FileExistsError:
            if pth.suffix == ".md5":
                # always reload md5 files
                pth.unlink()
                self.__download_file(ftp_, filename, tries + 1)
            elif not self.__checksum(pth):
                self.__download_file(ftp_, filename, tries + 1)
            else:
                logging.error('File %s already downloaded', filename)

    @staticmethod
    def get_names(dir_: Iterable) -> Set[str]:
        return set([Path(name).name for name in dir_])

    def __diff(
        self, source_dir: str, source_file_list: List[str], ext_: str
    ) -> Generator[Path, None, None]:
        pth = self.target_dir / source_dir
        targets = self.get_names(pth.glob(f"*.{ext_}"))
        sources = set(
            filter(
                lambda x: Path(x).suffix == f".{ext_}",
                self.get_names(source_file_list),
            )
        )
        diff_names = sources - targets

        if diff_names:
            ftp_dir = Path(source_file_list[0]).parent
            for name in diff_names:
                yield ftp_dir / name

    def full_download(self):
        source_dir = f"{self.target}/CURRENT-Full/{self.fmt}"

        with FTP(self.ftp, "anonymous", "anonymous@domain.com") as ftp_:
            # download md5
            for ext_ in self.sup_extensions:
                for remote_file in self.__diff(
                    source_dir, ftp_.list(source_dir), ext_[1:]
                ):
                    try:
                        logging.info("start download %s", remote_file)
                        self.__download_file(ftp_, remote_file)
                    except (ftplib.error_temp, EOFError, BaseException) as err_:
                        logging.warning("%s, sleeping 5 seconds", err_)
                        time.sleep(5)
                        logging.warning("retry...")
                        self.full_download()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Pubchem crawler')
    parser.add_argument('mode', type=str,
                        help='Modes: download')
    parser.add_argument('--pubchem-dir',
                        type=str,
                        default='./pubchem_dir',
                        help='Pubchem directory path')

    logging.basicConfig(level=logging.INFO)

    args = parser.parse_args()
    mode = args.mode
    pubchem_dir = Path(args.pubchem_dir)
    if mode == 'download':
        loader = PubchemLoader('compounds', pubchem_dir, 'sdf')
        loader.full_download()
    else:
        logging.error('Unknown mode %s', mode)
