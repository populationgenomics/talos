"""
Routes for analyst decisions — create, update, and view decision history.
"""

from fastapi import APIRouter, HTTPException, Request

from talos.web import database as db
from talos.web.models import DecisionHistoryEntry, DecisionRequest, DecisionResponse

router = APIRouter()


@router.post('/api/decisions', response_model=DecisionResponse)
async def api_upsert_decision(request: Request, body: DecisionRequest):
    conn = request.app.state.db
    try:
        db.upsert_decision(
            conn,
            sample=body.sample,
            chrom=body.chrom,
            pos=body.pos,
            ref=body.ref,
            alt=body.alt,
            gene=body.gene,
            status=body.status,
            acmg_class=body.acmg_class,
            analyst=body.analyst,
            comment=body.comment,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from None

    return db.get_decision(
        conn,
        sample=body.sample,
        chrom=body.chrom,
        pos=body.pos,
        ref=body.ref,
        alt=body.alt,
        gene=body.gene,
    )


@router.get('/api/decisions/{decision_id}/history', response_model=list[DecisionHistoryEntry])
async def api_decision_history(request: Request, decision_id: int):
    history = db.get_decision_history(request.app.state.db, decision_id)
    if not history:
        raise HTTPException(status_code=404, detail='Decision not found')
    return history
