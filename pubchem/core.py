import logging
from datetime import datetime

from ftpretty import ftpretty


def info(
    msg_: str,
) -> None:
    now_ = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{now_} [INFO] {msg_}")


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
