from enum import Enum


class DatasetSize(str, Enum):
    XSMALL = "xsmall"
    SMALL = "small"
    LARGE = "large"
    XLARGE = "xlarge"
