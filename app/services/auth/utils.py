import ipaddress
import os
import hmac
import hashlib
import secrets
import base64
from typing import Iterable, Optional, Tuple
from app.utils.logging_utils import logger


key_b64 = os.environ.get("APIKEY_HMAC_KEY_B64")
if key_b64:
    SERVER_HMAC_KEY = base64.b64decode(key_b64)
else:
    SERVER_HMAC_KEY = b""
    logger.warning("APIKEY_HMAC_KEY_B64 not set, HMAC operations will not work.")


def _b64url(nbytes: int) -> str:
    return base64.urlsafe_b64encode(secrets.token_bytes(nbytes)).rstrip(b"=").decode()


def generate_api_key(prefix_len_bytes: int = 5) -> Tuple[str, str, str]:
    """
    Returns (full_key, prefix).
    Format: sk_{env}_{prefixid}.{secret}
    secret: 32 bytes (~256-bit) base64url, no padding.
    """

    prefix_id = _b64url(prefix_len_bytes)  # short, non-secret lookup id
    secret = _b64url(32)
    full_key = f"{prefix_id}.{secret}"
    return full_key, prefix_id


def hmac_hash(full_key: str) -> bytes:
    digest = hmac.new(SERVER_HMAC_KEY, full_key.encode("utf-8"), hashlib.sha256).digest()
    return base64.b64encode(digest).decode("ascii")


def redact(full_key_or_secret: str) -> str:
    # Displays prefix + last4 only
    if "." in full_key_or_secret:
        prefix, secret = full_key_or_secret.split(".", 1)
        return f"{prefix}.…{secret[-4:]}"
    return "…" + full_key_or_secret[-4:]


def ip_allowed(client_ip: str, allowlist: Optional[Iterable[str]]) -> bool:
    if not allowlist:
        return True
    ip = ipaddress.ip_address(client_ip)
    return any(ip in ipaddress.ip_network(cidr) for cidr in allowlist)
