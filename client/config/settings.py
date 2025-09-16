from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """
    Application settings.
    """

    # Base URL of the backend API
    API_BASE_URL: str = "http://localhost:8000"

    # Optional API key or token
    API_KEY: str | None = None

    # Request timeout in seconds
    REQUEST_TIMEOUT: int = 30

    # Logging level
    LOG_LEVEL: str = "INFO"

    # class Config:
    #     env_file = ".env"
    #     env_file_encoding = "utf-8"


# singleton instance
settings = Settings()
