from app import models
from app.setup.database import SessionLocal


class AdminUser:
    def __init__(self):
        self._admin_user_id = None

    @property
    def admin_user_id(self):
        if self._admin_user_id is None:
            db = SessionLocal()
            try:
                admin = db.query(models.User).filter(models.User.email == "admin@datagrok.ai").first()
                if not admin:
                    raise Exception("Admin user not found.")
                self._admin_user_id = admin.id
            finally:
                db.close()
        return self._admin_user_id


admin = AdminUser()
