import jwt
from datetime import datetime, timedelta
from app.setup.config import SECRET_KEY, ALGORITHM, JWT_EXP_SECONDS


def create_jwt_token(user_id: int, session_id: int, role: str):
    now = datetime.now()
    payload = {
        "sub": str(user_id),
        "session_id": str(session_id),
        "role": role,
        "exp": now + timedelta(seconds=JWT_EXP_SECONDS),
    }
    token = jwt.encode(payload, SECRET_KEY, algorithm=ALGORITHM)
    return token


def decode_jwt_token(token: str):
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        return payload
    except jwt.JWTError:
        raise Exception("Invalid token")
