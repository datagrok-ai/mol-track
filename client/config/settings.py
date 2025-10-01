from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """
    Application settings.
    """

    API_BASE_URL: str = "http://localhost:8000"
    API_KEY: str | None = None
    REQUEST_TIMEOUT: int = 30

    # class Config:
    #     env_file = ".env"
    #     env_file_encoding = "utf-8"


# singleton instance
settings = Settings()
