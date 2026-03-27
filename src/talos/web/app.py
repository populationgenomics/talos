"""
Talos web application — FastAPI app factory and CLI entry point.
"""

import argparse
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from talos.web.database import init_db
from talos.web.services.result_loader import ResultCache

WEB_DIR = Path(__file__).parent
STATIC_DIR = WEB_DIR / 'static'


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialise database and result cache on startup."""
    app.state.db = init_db(app.state.db_path)
    app.state.result_cache = ResultCache()
    yield
    app.state.db.close()


def create_app(db_path: str = 'talos_decisions.db') -> FastAPI:
    app = FastAPI(title='Talos', lifespan=lifespan)
    app.state.db_path = db_path

    # include API route modules
    from talos.web.routes.decisions import router as decisions_router
    from talos.web.routes.runs import router as runs_router

    app.include_router(runs_router)
    app.include_router(decisions_router)

    # mount static files (CSS, JS, HTML pages) — must come after routes
    app.mount('/static', StaticFiles(directory=str(STATIC_DIR)), name='static')

    # serve static HTML pages for the SPA-style navigation
    @app.get('/')
    async def index() -> FileResponse:
        return FileResponse(STATIC_DIR / 'index.html')

    @app.get('/runs/{run_id}')
    async def run_detail(run_id: int) -> FileResponse:  # noqa: ARG001
        return FileResponse(STATIC_DIR / 'run.html')

    return app


def cli_main():
    parser = argparse.ArgumentParser(description='Talos web interface')
    parser.add_argument('--db-path', default='talos_decisions.db', help='Path to SQLite database file')
    parser.add_argument('--host', default='127.0.0.1', help='Host to bind to')
    parser.add_argument('--port', type=int, default=8080, help='Port to bind to')
    args = parser.parse_args()

    import uvicorn

    app = create_app(db_path=args.db_path)
    uvicorn.run(app, host=args.host, port=args.port)
