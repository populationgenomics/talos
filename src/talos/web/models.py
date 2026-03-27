"""
Pydantic models for the Talos web API requests and responses.
"""

from pydantic import BaseModel


class RegisterRunRequest(BaseModel):
    label: str
    json_path: str
    dataset: str = ''
    notes: str = ''


class RunResponse(BaseModel):
    id: int
    label: str
    dataset: str
    json_path: str
    run_date: str
    registered: str
    notes: str


class DecisionRequest(BaseModel):
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    status: str = 'pending'
    acmg_class: str | None = None
    analyst: str = ''
    comment: str = ''


class DecisionResponse(BaseModel):
    id: int
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    status: str
    acmg_class: str | None = None
    analyst: str = ''
    comment: str = ''
    created_at: str = ''
    updated_at: str = ''


class DecisionHistoryEntry(BaseModel):
    id: int
    decision_id: int
    status: str
    acmg_class: str | None = None
    analyst: str = ''
    comment: str = ''
    created_at: str = ''
