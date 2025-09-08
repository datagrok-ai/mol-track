import uuid
from sqlalchemy.orm import Session
from app.models import User, Session as SessionModel


def get_user_by_email(db: Session, user_email: str):
    return db.query(User).filter(User.email == user_email).first()


def create_session(db: Session, user_id: int):
    db_session = SessionModel(user_id=user_id)
    db.add(db_session)
    db.commit()
    db.refresh(db_session)
    return db_session


def get_session_by_id(db: Session, session_id: uuid.UUID):
    return db.query(SessionModel).filter(SessionModel.id == session_id).first()
