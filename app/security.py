import bcrypt
import secrets
import string


def generate_api_key():
    """Generate (key_id, secret) for a new key."""
    key_id = secrets.token_hex(8)  # public ID
    secret = "".join(secrets.choice(string.ascii_letters + string.digits) for _ in range(32))
    return key_id, secret


def hash_secret(secret: str) -> str:
    return bcrypt.hashpw(secret.encode(), bcrypt.gensalt()).decode()


def verify_secret(secret: str, secret_hash: str) -> bool:
    return bcrypt.checkpw(secret.encode(), secret_hash.encode())
