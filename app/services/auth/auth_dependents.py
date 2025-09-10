from datetime import datetime, timezone
import hmac
from sqlalchemy.orm import Session
from fastapi import Depends, HTTPException, Header, Request
from fastapi import status
from app.setup.database import get_db

from app.services.auth.utils import hmac_hash, ip_allowed
from app.services.auth.api_key_service import get_key_record
from app.utils import enums


# def enforce_https(request: Request):
#     if request.url.scheme != "https":
#         raise HTTPException(status_code=403, detail="HTTP is not allowed. Only HTTPS is accepted.")


def verify_api_key(
    request: Request,
    x_api_key: str = Header(..., alias="X-API-Key"),
    db: Session = Depends(get_db),
):
    if not x_api_key:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, detail="Missing API key")

    try:
        prefix, secret = x_api_key.split(".", 1)
    except ValueError:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, detail="Malformed API key")

    rec = get_key_record(db, prefix)
    if not rec or rec.status != enums.APIKeyStatus.active:
        raise HTTPException(
            status.HTTP_401_UNAUTHORIZED, detail="This API key does not exist or it isn't active anymore"
        )

    calc = hmac_hash(x_api_key)
    if not hmac.compare_digest(calc, rec.secret_hash):
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, detail="Invalid API key")

    expires_at = rec.expires_at
    if expires_at:
        now = datetime.now(timezone.utc)
        if expires_at < now:
            raise HTTPException(status_code=401, detail="API key expired")

    client_ip = request.client.host if request.client else "0.0.0.0"
    if not ip_allowed(client_ip, rec.ip_allowlist):
        raise HTTPException(status.HTTP_403_FORBIDDEN, detail="IP not allowed")

    return {"privileges": rec.privileges}


def require_privileges(*required: str):
    def dep(auth=Depends(verify_api_key)):
        if not set(required).intersection(auth["privileges"]):
            raise HTTPException(status.HTTP_403_FORBIDDEN, detail="You do not have privileges to access this resource")
        return auth

    return dep
