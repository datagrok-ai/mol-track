from datetime import datetime
from functools import lru_cache
from typing import Any, Dict, List, Optional
import uuid
from fastapi import status

from fastapi import HTTPException

from app.models import ApiKey
from sqlalchemy.orm import Session
from app.services.auth.utils import generate_api_key, hmac_hash
from app.utils import enums


def create_key(owner_id: uuid.UUID, privileges: List[enums.AuthPrivileges], ip_allowlist: List[str], db: Session):
    full_api_key, prefix = generate_api_key()
    rec = ApiKey(
        owner_id=owner_id,
        prefix=prefix,
        privileges=privileges,
        status=enums.APIKeyStatus.active.value,
        secret_hash=hmac_hash(full_api_key),
        created_at=datetime.now(),
        ip_allowlist=ip_allowlist,
    )
    try:
        db.add(rec)
        db.commit()
    except Exception as e:
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"An unexpected error occurred: {str(e)}"
        )

    return full_api_key


@lru_cache(maxsize=128)
def get_key_record(db: Session, prefix: str) -> Optional[Dict[str, Any]]:
    """
    Fetch key record by prefix.
    Returns dict with: stored_hash (str), status (str), scopes (list[str]),
    owner_id (UUID), ip_allowlist (list[str])
    """

    return db.query(ApiKey).filter(ApiKey.prefix == prefix).first()
