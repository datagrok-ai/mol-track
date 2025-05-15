import logging
from sqlalchemy import text
from sqlmodel import SQLModel
from sqlalchemy.orm import Session
from uuid import uuid4
from datetime import datetime
from .session import engine
from models import User

admin_user_id: str | None = None


def create_db_and_tables():
    SQLModel.metadata.create_all(engine)


def create_first_user_if_needed(session: Session):
    if session.query(User).first():
        return

    user = User(
        id=uuid4(),
        email="admin@example.com",
        first_name="Admin",
        last_name="User",
        has_password=True,
        is_active=True,
        created_at=datetime.utcnow(),
        updated_at=datetime.utcnow(),
        created_by=None,
        updated_by=None,
    )
    session.add(user)
    session.commit()
    session.refresh(user)


def init_db():
    with engine.connect() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS moltrack"))
        connection.commit()

    create_db_and_tables()

    with Session(engine) as session:
        create_first_user_if_needed(session)

        admin = session.query(User).filter(User.first_name == "Admin").first()
        if not admin:
            raise Exception("Admin user not found.")

        global admin_user_id
        admin_user_id = admin.id
        logging.info(f"Admin user ID: {admin_user_id}")
