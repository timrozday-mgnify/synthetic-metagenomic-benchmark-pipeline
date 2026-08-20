"""Interactive plots for reusable superresolution mis-mapping matrices."""

from __future__ import annotations

import math
from html import escape
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import HTML


_PALETTE = (
    "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
    "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
)
_SIZE = 820
_RADIUS = 322
_BAND = 13
_PAD_ANGLE = 0.006


def load_mismapping_matrix(path: Path, max_references: int = 80) -> tuple[pd.DataFrame, str, int]:
    """Load a dense or sparse matrix, retaining its most mis-mapped references.

    Sparse shotgun matrices use source/destination/probability triplets.  Limiting
    those to their largest off-diagonal masses avoids materialising a potentially
    very large chunk-by-chunk square frame merely to draw a report figure.
    """
    header = pd.read_csv(path, nrows=0).columns.tolist()
    sparse_columns = next(
        (columns for columns in (("src_chunk", "dst_chunk", "prob"), ("src", "dst", "prob"))
         if set(columns).issubset(header)),
        None,
    )
    if sparse_columns is None:
        matrix = pd.read_csv(path, index_col=0)
        matrix.index = matrix.index.astype(str)
        matrix.columns = matrix.columns.astype(str)
        matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        references = _select_references(matrix, max_references)
        return matrix.reindex(index=references, columns=references, fill_value=0.0), "amplicon", len(matrix)

    source_column, target_column, probability_column = sparse_columns
    triplets = pd.read_csv(
        path,
        usecols=[source_column, target_column, probability_column],
        dtype={source_column: str, target_column: str},
    ).rename(columns={source_column: "source", target_column: "target", probability_column: "probability"})
    triplets["probability"] = pd.to_numeric(triplets["probability"], errors="coerce").fillna(0.0)
    references = sorted(set(triplets["source"]) | set(triplets["target"]))
    retained = _select_sparse_references(triplets, references, max_references)
    matrix = (
        triplets[triplets["source"].isin(retained) & triplets["target"].isin(retained)]
        .pivot_table(index="source", columns="target", values="probability", aggfunc="sum")
        .reindex(index=retained, columns=retained, fill_value=0.0)
        .fillna(0.0)
    )
    return matrix, "shotgun", len(references)


def mismapping_heatmap(matrix: pd.DataFrame, title: str) -> go.Figure | None:
    """Return a log-scaled off-diagonal mis-mapping heatmap."""
    if matrix.empty:
        return None
    references = list(matrix.index.astype(str))
    values = matrix.to_numpy(dtype=float)
    masked = values.copy()
    np.fill_diagonal(masked, np.nan)
    with np.errstate(divide="ignore"):
        logged = np.log10(np.where(masked > 0, masked, np.nan))

    customdata = np.empty(values.shape + (5,), dtype=object)
    for i, source in enumerate(references):
        for j, target in enumerate(references):
            customdata[i, j] = (source, target, masked[i, j], masked[j, i], values[i, i])
    figure = go.Figure(go.Heatmap(
        z=logged,
        x=references,
        y=references,
        customdata=customdata,
        colorscale="Magma",
        zmin=-3,
        zmax=0,
        colorbar=dict(title="M[a,j]", tickvals=[-3, -2, -1, 0],
                      ticktext=["0.1%", "1%", "10%", "100%"]),
        hovertemplate=("From %{customdata[0]}<br>To %{customdata[1]}<br>"
                       "Mis-mapped: %{customdata[2]:.2%}<br>"
                       "Reverse: %{customdata[3]:.2%}<br>"
                       "Source keeps: %{customdata[4]:.1%}<extra></extra>"),
    ))
    figure.update_layout(
        title=title,
        xaxis=dict(title="Mapped to", showticklabels=False, constrain="domain"),
        yaxis=dict(title="Simulated from", showticklabels=False, autorange="reversed",
                   scaleanchor="x", constrain="domain"),
        template="simple_white",
        height=680,
    )
    return figure


def render_mismapping_chord(
    matrix: pd.DataFrame,
    title: str,
    min_value: float = 0.005,
) -> HTML | None:
    """Return the existing-report-style SVG ribbon plot for a mis-mapping matrix."""
    if matrix.empty:
        return None
    references = list(matrix.index.astype(str))
    values = matrix.to_numpy(dtype=float)
    flows = [
        (source, target, float(values[i, j]))
        for i, source in enumerate(references)
        for j, target in enumerate(references)
        if i != j and np.isfinite(values[i, j]) and values[i, j] >= min_value
    ]
    if not flows:
        return None
    flows.sort(key=lambda flow: flow[2])
    arcs = _chord_arcs(references, flows)
    colours = _reference_colours(references)
    spans = {
        reference: (arc["end"] - arc["start"])
        / sum(value for source, target, value in flows if reference in (source, target))
        for reference, arc in arcs.items()
    }
    ribbons = []
    for source, target, value in sorted(flows, key=lambda flow: (flow[0], flow[1])):
        source_span = _take(arcs[source], value * spans[source])
        target_span = _take(arcs[target], value * spans[target])
        ribbons.append((source, target, value, _ribbon_path(source_span, target_span)))

    index = {reference: i for i, reference in enumerate(references)}
    outgoing = {reference: 0.0 for reference in arcs}
    incoming = {reference: 0.0 for reference in arcs}
    for source, target, value in flows:
        outgoing[source] += value
        incoming[target] += value
    ribbon_svg = "\n".join(
        (f'<path class="ribbon" d="{path}" fill="{colours[_taxon(source)]}" '
         f'data-source="{escape(source)}" data-target="{escape(target)}" '
         f'data-value="{value}" data-reverse="{values[index[target], index[source]]}" '
         f'data-diagonal="{values[index[source], index[source]]}"></path>')
        for source, target, value, path in ribbons
    )
    arc_svg = "\n".join(
        (f'<path class="arc" d="{_arc_path(arc["start"], arc["end"])}" '
         f'fill="{colours[_taxon(reference)]}" data-reference="{escape(reference)}" '
         f'data-outgoing="{outgoing[reference]}" data-incoming="{incoming[reference]}" '
         f'data-diagonal="{values[index[reference], index[reference]]}"></path>')
        for reference, arc in arcs.items()
    )
    widget_id = f"mismapping-chord-{uuid4().hex}"
    legend = "".join(
        f'<span><i style="background:{colour}"></i>{escape(taxon)}</span>'
        for taxon, colour in colours.items()
    )
    return HTML(f"""
<style>
#{widget_id} {{ max-width:{_SIZE}px; position:relative; font-family:system-ui,sans-serif; }}
#{widget_id} svg {{ width:100%; height:auto; overflow:visible; }}
#{widget_id} .ribbon {{ cursor:pointer; opacity:.52; transition:opacity .1s; }}
#{widget_id} svg:hover .ribbon {{ opacity:.06; }}
#{widget_id} .ribbon:hover, #{widget_id} .ribbon.is-active {{ opacity:.95; stroke:#111; stroke-width:.6; }}
#{widget_id} .arc {{ cursor:pointer; stroke:#fff; stroke-width:.5; }}
#{widget_id} .tooltip {{ display:none; position:absolute; z-index:2; padding:8px 10px; background:#fffffff8; border:1px solid #ccc; border-radius:4px; font-size:12px; pointer-events:none; }}
#{widget_id} .legend {{ display:flex; flex-wrap:wrap; gap:4px 14px; font-size:12px; }}
#{widget_id} .legend i {{ display:inline-block; width:11px; height:11px; margin-right:4px; border-radius:2px; }}
</style>
<div id="{widget_id}"><h4>{escape(title)}</h4><svg viewBox="0 0 {_SIZE} {_SIZE}" role="img" aria-label="Mis-mapping between references">{ribbon_svg}{arc_svg}</svg><p>{len(arcs)} references; {len(ribbons)} off-diagonal mappings at least {min_value:.1%}. Hover a ribbon or arc for values.</p><div class="legend">{legend}</div><div class="tooltip"></div></div>
<script>(() => {{
const root=document.getElementById("{widget_id}"), tip=root.querySelector('.tooltip'), ribbons=[...root.querySelectorAll('.ribbon')];
const pct=value => {{ const n=Number(value)*100; return (n>=1?n.toFixed(1):n.toFixed(2))+'%'; }};
const move=event => {{ const b=root.getBoundingClientRect(); tip.style.left=(event.clientX-b.left+14)+'px'; tip.style.top=(event.clientY-b.top+14)+'px'; }};
ribbons.forEach(r => r.addEventListener('mousemove', event => {{ const d=r.dataset; tip.innerHTML='<b>'+d.source+' → '+d.target+'</b><br>Mis-mapped: '+pct(d.value)+'<br>Reverse: '+pct(d.reverse)+'<br>Source keeps: '+pct(d.diagonal); tip.style.display='block'; move(event); }}));
ribbons.forEach(r => r.addEventListener('mouseleave', () => tip.style.display='none'));
root.querySelectorAll('.arc').forEach(a => {{ a.addEventListener('mousemove', event => {{ const d=a.dataset; ribbons.forEach(r => r.classList.toggle('is-active', r.dataset.source===d.reference || r.dataset.target===d.reference)); tip.innerHTML='<b>'+d.reference+'</b><br>Leaks out: '+pct(d.outgoing)+'<br>Attracts in: '+pct(d.incoming)+'<br>Keeps its own: '+pct(d.diagonal); tip.style.display='block'; move(event); }}); a.addEventListener('mouseleave', () => {{ ribbons.forEach(r => r.classList.remove('is-active')); tip.style.display='none'; }}); }});
}})();</script>""")


def _select_references(matrix: pd.DataFrame, limit: int) -> list[str]:
    references = list(matrix.index.astype(str))
    if len(references) <= limit:
        return references
    values = matrix.reindex(index=references, columns=references).to_numpy(dtype=float)
    np.fill_diagonal(values, 0.0)
    masses = pd.Series(values.sum(axis=0) + values.sum(axis=1), index=references)
    return list(masses.nlargest(limit).index)


def _select_sparse_references(
    triplets: pd.DataFrame,
    references: list[str],
    limit: int,
) -> list[str]:
    if len(references) <= limit:
        return references
    off_diagonal = triplets[triplets["source"] != triplets["target"]]
    outgoing = off_diagonal.groupby("source")["probability"].sum()
    incoming = off_diagonal.groupby("target")["probability"].sum()
    masses = outgoing.add(incoming, fill_value=0.0).reindex(references, fill_value=0.0)
    return list(masses.nlargest(limit).index)


def _taxon(reference: str) -> str:
    return reference.split("|", 1)[0]


def _reference_colours(references: list[str]) -> dict[str, str]:
    taxa = sorted({_taxon(reference) for reference in references})
    return {taxon: _PALETTE[index % len(_PALETTE)] for index, taxon in enumerate(taxa)}


def _chord_arcs(
    references: list[str],
    flows: list[tuple[str, str, float]],
) -> dict[str, dict[str, float]]:
    outgoing: dict[str, float] = {}
    incoming: dict[str, float] = {}
    for source, target, value in flows:
        outgoing[source] = outgoing.get(source, 0.0) + value
        incoming[target] = incoming.get(target, 0.0) + value
    active = [reference for reference in references if reference in outgoing or reference in incoming]
    total = sum(outgoing.values()) + sum(incoming.values())
    available = 2 * math.pi - _PAD_ANGLE * len(active)
    angle = -math.pi / 2
    arcs: dict[str, dict[str, float]] = {}
    for reference in active:
        span = available * (outgoing.get(reference, 0.0) + incoming.get(reference, 0.0)) / total
        arcs[reference] = {"start": angle, "end": angle + span, "cursor": angle}
        angle += span + _PAD_ANGLE
    return arcs


def _take(arc: dict[str, float], span: float) -> tuple[float, float]:
    start = arc["cursor"]
    arc["cursor"] = start + span
    return start, arc["cursor"]


def _point(angle: float, radius: float) -> tuple[float, float]:
    centre = _SIZE / 2
    return centre + radius * math.cos(angle), centre + radius * math.sin(angle)


def _arc_path(start: float, end: float) -> str:
    outer = _RADIUS + _BAND
    large = int(end - start > math.pi)
    x0, y0 = _point(start, _RADIUS)
    x1, y1 = _point(end, _RADIUS)
    x2, y2 = _point(end, outer)
    x3, y3 = _point(start, outer)
    return (f"M{x0:.1f},{y0:.1f} A{_RADIUS},{_RADIUS} 0 {large} 1 {x1:.1f},{y1:.1f} "
            f"L{x2:.1f},{y2:.1f} A{outer},{outer} 0 {large} 0 {x3:.1f},{y3:.1f} Z")


def _ribbon_path(source_span: tuple[float, float], target_span: tuple[float, float]) -> str:
    centre = _SIZE / 2
    source_start, source_end = source_span
    target_start, target_end = target_span
    sx0, sy0 = _point(source_start, _RADIUS)
    sx1, sy1 = _point(source_end, _RADIUS)
    tx0, ty0 = _point(target_start, _RADIUS)
    tx1, ty1 = _point(target_end, _RADIUS)
    return (f"M{sx0:.1f},{sy0:.1f} A{_RADIUS},{_RADIUS} 0 0 1 {sx1:.1f},{sy1:.1f} "
            f"Q{centre},{centre} {tx0:.1f},{ty0:.1f} A{_RADIUS},{_RADIUS} 0 0 1 {tx1:.1f},{ty1:.1f} "
            f"Q{centre},{centre} {sx0:.1f},{sy0:.1f} Z")
