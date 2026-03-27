"""
API routes for run management and variant data.
"""

from fastapi import APIRouter, HTTPException, Request

from talos.web import database as db
from talos.web.models import RegisterRunRequest, RunResponse
from talos.web.services.prepare_data import prepare_metadata, prepare_run_context
from talos.web.services.variant_join import collect_variant_keys, join_variants_with_decisions

router = APIRouter()


@router.get('/api/runs', response_model=list[RunResponse])
async def api_list_runs(request: Request):
    runs = db.list_runs(request.app.state.db)
    return [dict(r) for r in runs]


@router.post('/api/runs', response_model=RunResponse)
async def api_register_run(request: Request, body: RegisterRunRequest):
    conn = request.app.state.db
    cache = request.app.state.result_cache

    try:
        result_data = cache.get(body.json_path)
    except FileNotFoundError:
        raise HTTPException(status_code=400, detail=f'Result file not found: {body.json_path}') from None

    run_date = result_data.metadata.run_datetime

    try:
        run_id = db.register_run(
            conn,
            label=body.label,
            json_path=body.json_path,
            dataset=body.dataset,
            run_date=run_date,
            notes=body.notes,
        )
    except Exception as e:
        if 'UNIQUE constraint' in str(e):
            raise HTTPException(status_code=409, detail='This JSON file is already registered') from None
        raise

    run = db.get_run(conn, run_id)
    return dict(run)


@router.delete('/api/runs/{run_id}')
async def api_delete_run(request: Request, run_id: int):
    if not db.delete_run(request.app.state.db, run_id):
        raise HTTPException(status_code=404, detail='Run not found')
    return {'ok': True}


@router.get('/api/runs/{run_id}/data')
async def api_run_data(request: Request, run_id: int):
    """Return full run data: metadata, samples, and variants with decisions."""
    conn = request.app.state.db
    run = db.get_run(conn, run_id)
    if not run:
        raise HTTPException(status_code=404, detail='Run not found')

    cache = request.app.state.result_cache
    try:
        result_data = cache.get(run['json_path'])
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f'Result file not found: {run["json_path"]}') from None

    keys = collect_variant_keys(result_data)
    decisions = db.get_decisions_for_keys(conn, keys)
    joined = join_variants_with_decisions(result_data, decisions)
    samples = prepare_run_context(result_data, joined)
    metadata = prepare_metadata(result_data)

    return {
        'run': dict(run),
        'metadata': metadata,
        'samples': samples,
    }
