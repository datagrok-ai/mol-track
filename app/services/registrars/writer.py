import io
import csv
import orjson


class StreamingResultWriter:
    def __init__(self, format: str):
        self.format = format
        self._header_written = False
        self._first_json_row = True

    def stream_rows(self, rows_iter):
        if self.format == "csv":
            yield from self._stream_csv(rows_iter)
        if self.format == "json":
            yield from self._stream_json(rows_iter)
        else:
            yield from self._stream_sdf(rows_iter)

    def _stream_csv(self, rows_iter):
        output = io.StringIO()
        writer = None

        for rows in rows_iter:
            if not writer:
                writer = csv.DictWriter(output, fieldnames=rows[0].keys())
                writer.writeheader()
                yield output.getvalue()
                output.seek(0)
                output.truncate(0)

            writer.writerows(rows)
            yield output.getvalue()
            output.seek(0)
            output.truncate(0)

    def _stream_json(self, rows_iter):
        yield "["
        first = True
        for rows in rows_iter:
            for row in rows:
                if not first:
                    yield ","
                else:
                    first = False
                yield orjson.dumps(row).decode("utf-8")
        yield "]"

    def _stream_sdf(self, rows_iter):
        for rows in rows_iter:
            for row in rows:
                molfile = row.get("original_molfile", "")
                yield molfile + "\n"
                for key, value in row.items():
                    if key == "original_molfile":
                        continue
                    yield f"> <{key}>\n{value}\n\n"
                yield "$$$$\n"
