import os
import csv
import orjson
from app.utils import enums


class ResultWriter:
    def __init__(self, output_path: str, format: str = enums.OutputFormat.json):
        if not output_path:
            raise ValueError("Output path must be specified")

        self.output_path = output_path
        self.format = format

        if os.path.exists(self.output_path):
            os.remove(self.output_path)

        self._file = open(self.output_path, "a", newline="", encoding="utf-8")
        self._writer = None
        self._header_written = False
        self._first_json_row = True

        if self.format == enums.OutputFormat.json:
            self._file.write("[")

    def write_rows(self, rows):
        if not rows:
            return

        if self.format == enums.OutputFormat.csv:
            if self._writer is None:
                self._writer = csv.DictWriter(self._file, fieldnames=rows[0].keys())

            if not self._header_written:
                self._writer.writeheader()
                self._header_written = True

            self._writer.writerows(rows)

        elif self.format == enums.OutputFormat.json:
            for row in rows:
                if not self._first_json_row:
                    self._file.write(",")
                else:
                    self._first_json_row = False
                self._file.write(orjson.dumps(row).decode("utf-8"))

        self._file.flush()

    def close(self):
        if self._file:
            if self.format == enums.OutputFormat.json:
                self._file.write("]")

            self._file.close()
            self._file = None
            self._writer = None

    def delete(self):
        self.close()
        if self.output_path and os.path.exists(self.output_path):
            os.remove(self.output_path)
