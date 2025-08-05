import os
import csv


class ResultWriter:
    def __init__(self, output_path: str):
        if not output_path:
            raise ValueError("Output path must be specified")

        self.output_path = output_path

        if os.path.exists(self.output_path):
            os.remove(self.output_path)

        self._file = open(self.output_path, "a", newline="", encoding="utf-8")
        self._writer = None
        self._header_written = False

    def write_rows(self, rows):
        if not rows:
            return

        if self._writer is None:
            self._writer = csv.DictWriter(self._file, fieldnames=rows[0].keys())

        if not self._header_written:
            self._writer.writeheader()
            self._header_written = True

        self._writer.writerows(rows)
        self._file.flush()

    def close(self):
        if self._file:
            self._file.close()
            self._file = None
            self._writer = None

    def delete(self):
        self.close()
        if self.output_path and os.path.exists(self.output_path):
            os.remove(self.output_path)
