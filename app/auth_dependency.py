from fastapi import Depends, Header, HTTPException
from sqlalchemy.orm import Session
from app.models import ApiKey
from app.security import verify_secret
import datetime


def require_auth_scopes(required: list[str]):
    from app.main import get_db

    def checker(
        x_api_key: str = Header(..., alias="X-API-Key"),
        db: Session = Depends(get_db),
    ):
        try:
            key_id, secret = x_api_key.split(".", 1)
        except ValueError:
            raise HTTPException(status_code=401, detail="Malformed API key")

        api_key: ApiKey = db.query(ApiKey).filter(ApiKey.key_id == key_id).first()
        if not api_key or api_key.revoked:
            raise HTTPException(status_code=401, detail="Invalid or revoked key")

        if api_key.expires_at and api_key.expires_at < datetime.utcnow():
            raise HTTPException(status_code=401, detail="Key expired")

        if not verify_secret(secret, api_key.secret_hash):
            raise HTTPException(status_code=401, detail="Invalid key")

        scopes = api_key.scopes.split(",") if api_key.scopes else []
        if not any(s in scopes for s in required):
            raise HTTPException(status_code=403, detail="Not enough privileges")

        return scopes

    return checker
