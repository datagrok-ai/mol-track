from sqlalchemy.orm import Session
from app.setup.database import SessionLocal
from app.models import ApiKey
from app.security import generate_api_key, hash_secret
from app.utils import enums


def create_key(name: str, scopes: list[str]):
    db: Session = SessionLocal()
    key_id, secret = generate_api_key()
    secret_hash = hash_secret(secret)

    api_key = ApiKey(key_id=key_id, secret_hash=secret_hash, scopes=",".join(scopes))
    db.add(api_key)
    db.commit()
    db.refresh(api_key)

    print("== API KEY GENERATED ==")
    print(f"API key for {name}:\nX-API-Key: {key_id}.{secret}")
    print(f"Scopes: {scopes}")


if __name__ == "__main__":
    create_key("reader", [enums.AuthScopes.READER])
    create_key("writer", [enums.AuthScopes.WRITER])
    create_key("admin", [enums.AuthScopes.ADMIN])
