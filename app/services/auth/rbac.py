from fastapi import Depends, HTTPException, status


def require_role(required_role: str):
    def role_dependency(user=Depends("app.auth.dependencies.get_current_user")):
        if user["role"] != required_role:
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Insufficient role")
        return user

    return role_dependency
