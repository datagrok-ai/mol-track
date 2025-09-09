import os
import hmac
import hashlib
import secrets
import base64
from typing import Tuple


SERVER_HMAC_KEY = base64.b64decode(os.environ["APIKEY_HMAC_KEY_B64"])


def _b64url(nbytes: int) -> str:
    return base64.urlsafe_b64encode(secrets.token_bytes(nbytes)).rstrip(b"=").decode()


def generate_api_key(prefix_len_bytes: int = 5) -> Tuple[str, str, str]:
    """
    Returns (full_key, prefix, last4).
    Format: sk_{env}_{prefixid}.{secret}
    secret: 32 bytes (~256-bit) base64url, no padding.
    """

    prefix_id = _b64url(prefix_len_bytes)  # short, non-secret lookup id
    secret = _b64url(32)
    full_key = f"{prefix_id}.{secret}"
    last4 = secret[-4:]
    return full_key, prefix_id, last4


def hmac_hash(full_key: str) -> bytes:
    digest = hmac.new(SERVER_HMAC_KEY, full_key.encode("utf-8"), hashlib.sha256).digest()
    return base64.b64encode(digest).decode("ascii")


def redact(full_key_or_secret: str) -> str:
    # Displays prefix + last4 only
    if "." in full_key_or_secret:
        prefix, secret = full_key_or_secret.split(".", 1)
        return f"{prefix}.…{secret[-4:]}"
    return "…" + full_key_or_secret[-4:]
