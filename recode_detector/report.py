"""Plain-English PDF report for the database comparison. Needs reportlab (pure Python)."""
from __future__ import annotations

import datetime as _dt
import io
from collections import Counter

from . import __version__

TEAL = "#1f3a93"
AMBER = "#d9452b"
BLUE = "#2f7d6d"
GREY = "#6d6657"
LIGHT = "#f8f2e4"


def _need_reportlab():
    try:
        import reportlab  # noqa: F401
    except ImportError as e:
        raise ValueError("PDF reports need the reportlab package: pip install reportlab "
                         "(the Mac app installs it itself the next time it starts).") from e


def _n(v: float) -> str:
    return f"{int(round(v)):,}"


def _fi(v: float) -> str:
    return f"{v / 1e6:.1f} million" if v >= 1e6 else f"{v / 1e3:.0f} thousand" if v >= 1e3 else f"{v:.0f}"


def _esc(t: str) -> str:
    return str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _verdict_sentence(c: dict, ra: str, aa: str) -> str:
    alt, ref, both, none = c["alternative"], c["reference"], c["both"], c["none"]
    judged = alt + ref + both
    if judged == 0:
        return ("No position could be judged: none of the changed positions was covered by an identified peptide, "
                "so this run neither supports nor argues against the alternative version.")
    if alt and not ref and not both:
        return (f"Where the run can speak at all ({alt} positions), it saw only the alternative version ({aa}) and "
                f"never the standard one ({ra}). That is consistent with the alternative reading, but it rests on few "
                "peptides from a single run and needs confirming.")
    if ref and not alt and not both:
        return (f"Where the run can speak at all ({ref} positions), it saw only the standard version ({ra}). "
                "That does not support the alternative reading at those positions.")
    return (f"The run gives mixed evidence: the alternative version was seen alone at {alt} positions, the standard "
            f"version alone at {ref}, and both at {both}. Each position needs to be looked at individually.")


def spectra_summary(spectra: dict | None, sites: list[dict]) -> dict | None:
    """How the raw-spectrum check relates to what the search software reported."""
    if not spectra or not spectra.get("results"):
        return None
    obs = {(i, v): {p["sequence"] for p in s[k] if not p["also_elsewhere"]}
           for i, s in enumerate(sites) for v, k in (("alternative", "alt_peptides"), ("reference", "ref_peptides"))}
    out = {"alt_obs": 0, "alt_supported": 0, "alt_weak": 0, "alt_not": 0, "ref_ctr_supported": 0, "ref_ctr_weak": 0, "n": len(spectra["results"])}
    for x in spectra["results"]:
        site, ver, seq = x.get("site"), x.get("version"), x.get("peptide")
        if site is None:
            continue
        observed = seq in obs.get((site, ver), set())
        if observed and ver == "alternative":
            out["alt_obs"] += 1
            out["alt_supported" if x["verdict"] == "supported" else "alt_weak" if x["verdict"] == "weak" else "alt_not"] += 1
        elif not observed and ver == "reference" and not obs.get((site, "reference")):
            if x["verdict"] == "supported":
                out["ref_ctr_supported"] += 1
            elif x["verdict"] == "weak":
                out["ref_ctr_weak"] += 1
    out["conflict"] = out["alt_not"] > 0 or out["ref_ctr_supported"] > 0
    return out


def build_pdf(res, files: dict, run_date: _dt.date | None = None, note: str | None = None,
              qc: dict | None = None, spectra: dict | None = None) -> bytes:
    """res: dbcompare.Result. files: {slot: filename}. note: optional banner under the title (for example data)."""
    _need_reportlab()
    from reportlab.graphics.charts.barcharts import VerticalBarChart
    from reportlab.graphics.shapes import Drawing, Line, Rect, String
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,
                                    TableStyle)

    cov, summ = res.coverage, res.summary
    c = summ["status"]
    sites = res.sites
    changes = Counter((s["ref_aa"], s["alt_aa"]) for s in sites)
    (ra, aa), _ = changes.most_common(1)[0]
    one_change = len(changes) == 1
    change_txt = f"{ra} to {aa}" if one_change else "several kinds of change"
    run_date = run_date or _dt.date.today()

    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["BodyText"], fontName="Helvetica", fontSize=10.3, leading=14.6,
                          alignment=TA_LEFT, spaceAfter=6)
    small = ParagraphStyle("small", parent=body, fontSize=8.6, leading=11.5, textColor=colors.HexColor(GREY))
    h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold", fontSize=19, leading=23,
                        textColor=colors.HexColor(TEAL), spaceBefore=4, spaceAfter=8, keepWithNext=1)
    h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold", fontSize=13, leading=17,
                        textColor=colors.HexColor("#0f172a"), spaceBefore=12, spaceAfter=5, keepWithNext=1)
    title = ParagraphStyle("title", parent=h1, fontSize=26, leading=30, spaceAfter=4)
    bullet = ParagraphStyle("bullet", parent=body, leftIndent=14, bulletIndent=2, spaceAfter=3)
    cell = ParagraphStyle("cell", parent=body, fontSize=8.6, leading=11, spaceAfter=0)
    cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")

    def P(t, st=body):
        return Paragraph(t, st)

    def bullets(items):
        return [Paragraph(t, bullet, bulletText="•") for t in items]

    def chart(values, labels, w=250, h=140, color=TEAL, ylab="", fmt=lambda v: f"{v:g}", bar_colors=None):
        d = Drawing(w, h + 36)
        bc = VerticalBarChart()
        bc.x, bc.y, bc.width, bc.height = 38, 26, w - 48, h - 20
        bc.data = [list(values)]
        bc.categoryAxis.categoryNames = labels
        bc.categoryAxis.labels.fontSize = 7.2
        bc.categoryAxis.labels.dy = -2
        bc.valueAxis.labels.fontSize = 7.2
        bc.valueAxis.valueMin = 0
        bc.valueAxis.valueMax = max(max(values) * 1.15, 1)
        bc.valueAxis.visibleGrid = True
        bc.valueAxis.gridStrokeColor = colors.HexColor("#e2e8f0")
        bc.bars[0].fillColor = colors.HexColor(color)
        bc.bars[0].strokeColor = None
        bc.categoryAxis.labels.fontName = bc.valueAxis.labels.fontName = bc.barLabels.fontName = "Helvetica"
        for i, col in enumerate(bar_colors or []):
            bc.bars[(0, i)].fillColor = colors.HexColor(col)
        bc.barLabelFormat = fmt
        bc.barLabels.fontSize = 6.8
        bc.barLabels.nudge = 6
        d.add(bc)
        if ylab:
            d.add(String(0, h + 22, ylab, fontSize=7.4, fontName="Helvetica", fillColor=colors.HexColor(GREY)))
        return d

    def captioned(drawing, text):
        return [drawing, P(text, small)]

    def table(rows, widths, header=True, zebra=True):
        t = Table(rows, colWidths=widths, repeatRows=1 if header else 0)
        st = [("FONTNAME", (0, 0), (-1, -1), "Helvetica"), ("FONTSIZE", (0, 0), (-1, -1), 8.6),
              ("VALIGN", (0, 0), (-1, -1), "TOP"), ("BOTTOMPADDING", (0, 0), (-1, -1), 4), ("TOPPADDING", (0, 0), (-1, -1), 4),
              ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#cbd5e1"))]
        if header:
            st += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(TEAL)), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                   ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold")]
        if zebra:
            st += [("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1), [colors.white, colors.HexColor("#f8fafc")])]
        t.setStyle(TableStyle(st))
        return t

    story = []
    # ------------------------------------------------------------------ title
    story += [P("Database comparison report", title),
              P(f"Which version of each protein did this mass spectrometry run see? ({_esc(change_txt)})", ParagraphStyle(
                  "sub", parent=body, fontSize=12.5, leading=16, textColor=colors.HexColor(GREY))),
              P(f"Prepared {run_date.strftime('%d %B %Y')} with Recode Detector {__version__}. "
                f"Run analysed: {_esc(', '.join(res.samples))}.", small)]
    if note:
        story += [P(f"<b>{_esc(note)}</b>", ParagraphStyle("note", parent=body, fontSize=9.6, textColor=colors.HexColor(AMBER)))]
    story += [Spacer(1, 6)]

    # ------------------------------------------------------------------ short version
    n_judged = c["alternative"] + c["reference"] + c["both"]
    ss_ = spectra_summary(spectra, sites)
    short = bullets([
        f"The run identified <b>{_n(cov['n_peptides'])}</b> different peptides (short pieces of protein) and found at least one "
        f"for <b>{_n(cov['n_detected'])} of the {_n(cov['n_proteins'])}</b> proteins in the reference database "
        f"(<b>{100 * cov['n_detected'] / cov['n_proteins']:.0f}%</b>).",
        f"Across all proteins, about <b>{cov['overall_cov']:.0f}%</b> of the letters (amino acids) were covered by an identified peptide. "
        f"A typical detected protein had about {cov['median_cov']:.0f}% of its sequence seen.",
        f"The alternative database changes <b>{summ['n_sites']}</b> positions in <b>{summ['n_proteins']}</b> proteins "
        f"({_esc(change_txt)}).",
        f"The run could judge <b>{n_judged} of {summ['n_sites']}</b> of those positions: <b>{c['alternative']}</b> showed only the "
        f"alternative version, <b>{c['reference']}</b> only the standard version, <b>{c['both']}</b> both. "
        f"The other <b>{c['none']}</b> were not covered by any identified peptide, so the run says nothing about them."]
        + ([f"<b>Raw-spectrum check:</b> of the {ss_['alt_obs']} alternative-version peptides the search software identified, the actual "
            f"fragmentation spectra confirm <b>{ss_['alt_supported']}</b> and partly support {ss_['alt_weak']}; {ss_['alt_not']} could not be "
            f"confirmed. At <b>{ss_['ref_ctr_supported']}</b> positions the spectra instead show convincing evidence for the "
            f"<i>reference</i> version (the original letter), which the search did not report."] if ss_ else []))
    box = Table([[ [P("The short version", h2)] + short + [Spacer(1, 4),
                    P("<b>Bottom line.</b> " + _verdict_sentence(c, ra, aa)
                      + (" <b>However, the raw spectra do not fully agree with the search results (see \"Looking at the raw spectra\"), so "
                         "treat the conclusion as open until the key peptides have been checked in a spectrum viewer.</b>"
                         if ss_ and ss_["conflict"] else ""), body)] ]], colWidths=[6.9 * inch])
    box.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(LIGHT)),
                             ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor(TEAL)),
                             ("LEFTPADDING", (0, 0), (-1, -1), 12), ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                             ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
    story += [box]

    # ------------------------------------------------------------------ what was measured
    story += [P("What was measured", h2),
              P("The sample's proteins were cut into short pieces called <b>peptides</b> (with an enzyme, usually trypsin, that "
                "cuts after the letters K and R). A mass spectrometer measured the pieces, and software (DIA-NN) matched "
                "them against a list of protein sequences, a <b>database</b>, to work out which proteins were present and "
                "how much of each. A peptide is only identified if its sequence is in the database that was searched."),
              P(f"Here the search used a database that held both the <b>standard</b> protein sequences and "
                f"<b>alternative</b> versions of {summ['n_proteins']} proteins, written as if the genetic code were read "
                f"differently ({_esc(change_txt)}). This lets the data \"vote\" for one version or the other.")]

    # ------------------------------------------------------------------ coverage
    story += [P("How much of the organism was seen?", h2),
              P(f"This is called <b>coverage</b>. Two kinds matter. <b>Protein coverage</b> asks whether a protein was seen at all: "
                f"{_n(cov['n_detected'])} of {_n(cov['n_proteins'])} ({100 * cov['n_detected'] / cov['n_proteins']:.0f}%). "
                f"For {_n(cov['n_ge2'])} proteins at least two different peptides were found, which makes the identification more "
                f"dependable, and for {_n(cov['n_ge5'])} at least five. <b>Sequence coverage</b> asks how much of each protein was seen. "
                f"Even a well-detected protein is usually only partly covered, because some peptides are too short, too long or "
                f"too hard to detect. Here the middle protein had {cov['median_cov']:.0f}% of its letters covered, and "
                f"{_n(cov['n_cov50'])} proteins had at least half.")]
    hist = chart(cov["hist"], [f"{i * 10}+" for i in range(10)], ylab="Number of proteins")
    bl = cov["by_length"]
    bylen = chart([100 * b["detected"] / b["n"] if b["n"] else 0 for b in bl], [f"{b['label']}" for b in bl],
                  color=BLUE, ylab="% of proteins detected", fmt=lambda v: f"{v:.0f}%")
    story += [Table([[captioned(hist, "Figure 1. For the proteins that were detected: how much of each protein's sequence was covered "
                                      "(each bar starts at the percent shown, so 30+ means 30 to under 40)."),
                      captioned(bylen, "Figure 2. Share of proteins detected, by protein size (number of amino acids). "
                                       "Short proteins give few peptides, so they are missed more often.")]],
                    colWidths=[3.45 * inch, 3.45 * inch])]
    story += [P(f"In plain terms: the run saw most of the proteins, but only about a quarter of the total sequence. "
                f"No single run sees everything, and that limit decides what can be said about any one position, as the next sections show. "
                f"Peptide signal strength spanned about {cov['intensity_orders']:.1f} orders of magnitude, a wide range, "
                f"so weak proteins can sit near the detection limit. Peptides had a median of "
                f"{_median_len(cov)} amino acids, and {100 * cov['missed'].get('0', 0) / max(cov['n_peptides'], 1):.0f}% had no missed "
                f"cut by the enzyme.")]

    # ------------------------------------------------------------------ the question
    story += [P("The question: which version of the protein?", h2),
              P("At a changed position, the standard database has one letter (here "
                f"<b>{_esc(ra)}</b>) and the alternative has another (<b>{_esc(aa)}</b>). Only a peptide that passes "
                "<i>through</i> that position carries the answer: it will contain either one letter or the other. Peptides "
                "elsewhere in the protein are identical in both versions and say nothing. So for every changed position the "
                "program looks for identified peptides that span it, and counts one only if it does not also occur "
                "somewhere else in the database."),
              P("A side effect worth knowing: trypsin cuts after R but not after W. If a position really is W, the peptide that "
                "would have ended there is instead joined to the next piece, so the alternative peptide is usually longer than the "
                "standard one. That is one reason a position can go uncovered.")]
    unc = [s_ for s_ in sites if s_["status"] == "none"]
    if unc:
        both_ok = sum(1 for s_ in unc if s_["expected_ref"]["detectable"] and s_["expected_alt"]["detectable"])
        ref_only = sum(1 for s_ in unc if s_["expected_ref"]["detectable"] and not s_["expected_alt"]["detectable"])
        alt_only = sum(1 for s_ in unc if s_["expected_alt"]["detectable"] and not s_["expected_ref"]["detectable"])
        neither = len(unc) - both_ok - ref_only - alt_only
        story += [P(f"<b>Why were {len(unc)} positions not covered?</b> For each one the program works out which peptides the enzyme "
                    f"would produce and whether a search could see them (between 7 and 30 amino acids). At <b>{both_ok}</b> of the uncovered "
                    f"positions <i>both</i> versions would have made a detectable peptide, so the run simply did not identify one. "
                    f"At <b>{ref_only}</b> only the standard version would be detectable (the alternative peptide would be too long or too short), "
                    f"at <b>{alt_only}</b> only the alternative version, and at <b>{neither}</b> neither. Where a version could never be "
                    f"detected, its absence is not evidence about the biology.")]
    bars = chart([c["alternative"], c["reference"], c["both"], c["none"]],
                 ["Alternative", "Standard", "Both", "No peptide"], w=300, ylab="Number of positions",
                 bar_colors=[AMBER, TEAL, "#c98a1b", "#a39b8a"])
    story += [KeepTogether(captioned(bars, f"Figure 3. What the run showed at the {summ['n_sites']} changed positions."))]

    story += [P("Positions where the run could judge", h2)]
    shown = [s for s in sites if s["status"] != "none"]
    if shown:
        rows = [[P("Protein", cellb), P("Position", cellb), P("Change", cellb), P("Verdict", cellb),
                 P("Peptide seen", cellb), P("Signal", cellb)]]
        for s in shown:
            ps = [p for k in ("alt_peptides", "ref_peptides") for p in s[k] if not p["also_elsewhere"]]
            pep = "<br/>".join(f"{_esc(p['sequence'])} (residues {p['start']}-{p['end']})" for p in ps[:3])
            rows.append([P(_esc(_short(s["ref_protein"])), cell), P(str(s["position"]), cell),
                         P(f"{_esc(s['ref_aa'])} to {_esc(s['alt_aa'])}", cell), P(_verdict_word(s["status"]), cell),
                         P(pep, cell), P(_fi(sum(p["total_intensity"] for p in ps)), cell)])
        story += [table(rows, [1.1 * inch, 0.65 * inch, 0.7 * inch, 1.15 * inch, 2.6 * inch, 0.7 * inch])]
        story += [P("\"Signal\" is the summed measured intensity of the peptides shown, in arbitrary instrument units; it "
                    "shows relative strength, not an amount.", small)]
    else:
        story += [P("No changed position was covered by an identified peptide.")]
    story += [P("Comparison with a search against the standard database", h2)]
    if summ.get("has_b"):
        sc = res.scatter
        story += [P(f"The same data were also searched against the standard database alone (no alternative proteins). "
                    f"<b>{_n(sc['n_shared'])}</b> peptides were found by both searches, <b>{_n(sc['n_only_a'])}</b> only in the main search "
                    f"(which included the alternative proteins) and <b>{_n(sc['n_only_b'])}</b> only in the standard-only search. "
                    f"For the peptides found by both, the signal strengths agree closely "
                    f"(correlation {sc['r']:.2f}), as they should for the same sample. That suggests the two searches saw the same "
                    f"material and that adding the alternative proteins left the rest of the results largely unchanged.")]
        d = Drawing(250, 200)
        lo = min(min(sc["a"]), min(sc["b"])) - 0.5 if sc["a"] else 0
        hi = max(max(sc["a"]), max(sc["b"])) + 0.5 if sc["a"] else 1
        x0, y0, w, h = 38, 30, 200, 150
        d.add(Rect(x0, y0, w, h, strokeColor=colors.HexColor("#94a3b8"), fillColor=None, strokeWidth=0.6))
        d.add(Line(x0, y0, x0 + w, y0 + h, strokeColor=colors.HexColor(AMBER), strokeDashArray=[3, 3], strokeWidth=0.8))
        from reportlab.graphics.shapes import Circle
        step = max(1, len(sc["a"]) // 1200)
        for xa, yb in list(zip(sc["a"], sc["b"]))[::step]:
            d.add(Circle(x0 + (xa - lo) / (hi - lo) * w, y0 + (yb - lo) / (hi - lo) * h, 1.3,
                         fillColor=colors.HexColor(TEAL), strokeColor=None, fillOpacity=0.45))
        for t in range(int(lo) + 1, int(hi), max(1, int((hi - lo) // 5))):
            d.add(String(x0 + (t - lo) / (hi - lo) * w - 5, y0 - 10, str(t), fontSize=7, fontName="Helvetica"))
            d.add(String(x0 - 16, y0 + (t - lo) / (hi - lo) * h - 2, str(t), fontSize=7, fontName="Helvetica"))
        d.add(String(x0 + 20, 6, "main search (with alternative proteins)", fontSize=7.4, fontName="Helvetica",
                     fillColor=colors.HexColor(GREY)))
        d.add(String(0, y0 + h + 8, "standard-database search", fontSize=7.4, fontName="Helvetica", fillColor=colors.HexColor(GREY)))
        story += [KeepTogether(captioned(d, "Figure 4. Each dot is a peptide found by both searches (log2 of its signal). "
                                            "Dots near the dashed line mean the two searches agree."))]
        hits = [s for s in sites if s["ref_peptides_b"]]
        if hits:
            story += [P(f"At <b>{len(hits)}</b> of the changed positions, the standard-only search identified a peptide carrying the "
                        f"original letter ({_esc(ra)}):")]
            rows = [[P("Protein", cellb), P("Position", cellb), P("Main-search verdict", cellb),
                     P("Standard-search peptide", cellb), P("Also in main search?", cellb)]]
            for s_ in hits:
                rows.append([P(_esc(_short(s_["ref_protein"])), cell), P(str(s_["position"]), cell),
                             P(_verdict_word(s_["status"]), cell),
                             P("<br/>".join(_esc(p["sequence"]) for p in s_["ref_peptides_b"][:3]), cell),
                             P("<br/>".join("yes" if p["seen_in_a"] else "no" for p in s_["ref_peptides_b"][:3]), cell)])
            story += [table(rows, [1.2 * inch, 0.7 * inch, 1.6 * inch, 2.3 * inch, 1.1 * inch]),
                      P("How to read this: the standard-only search has no alternative proteins, so it can only explain a spectrum with the "
                        "original letter. If a peptide with the original letter shows up there but is absent from the main search, "
                        "while the main search found the alternative version at the same position, the data prefer the alternative. "
                        "It can also just mean the two searches differ in what they accept, so look at the spectra before relying on it.",
                        small)]
        else:
            story += [P("The standard-only search did not identify a peptide spanning any changed position, so it adds no further "
                        "information about them.")]
    else:
        story += [P("<b>Not included in this report.</b> The data were searched once, against a database that held both the standard "
                    "and the alternative proteins, so there is nothing to compare it with yet."),
                  P("The comparison is worth adding. Searching the same raw data again against the <i>standard</i> database alone "
                    "(without the alternative proteins) shows what the data look like when only the original reading is allowed. "
                    "Peptides that carry the original letter at a changed position can only be found in that search, and a position "
                    "where they appear there but not in the main search may be one where the data favour the alternative. "
                    "To add it: run the same files through DIA-NN with the standard FASTA only, then load the new "
                    "<font face='Courier'>report.pr_matrix.tsv</font> in the program's fourth file box and download this report again.")]

    # ------------------------------------------------------------------ quality checks
    if qc and qc.get("checks"):
        story += [P("Quality checks on the run", h2),
                  P("These checks ask whether the run itself looks healthy, using rules of thumb that are shown next to each result. "
                    "A warning is a reason to look closer, not proof of a problem.")]
        word = {"pass": "Pass", "warn": "Check", "fail": "Problem", "info": "Info"}
        col = {"pass": "#2f7d6d", "warn": "#c98a1b", "fail": "#d9452b", "info": "#6d6657"}
        rows = [[P("Check", cellb), P("Result", cellb), P("Verdict", cellb), P("What it means", cellb)]]
        for c_ in qc["checks"]:
            rows.append([P(_esc(c_["name"]), cell), P(_esc(c_["value"]), cell),
                         P(f"<font color='{col[c_['status']]}'><b>{word[c_['status']]}</b></font> ({_esc(c_['rule'])})", cell),
                         P(_esc(c_["why"]), cell)])
        story += [table(rows, [1.9 * inch, 1.1 * inch, 1.5 * inch, 2.4 * inch])]

    # ------------------------------------------------------------------ spectrum evidence
    if spectra and spectra.get("results"):
        story += _spectra_section(spectra, sites, P, h2, cell, cellb, small, table, captioned, colors, inch, Drawing, Line, String, Rect, KeepTogether)

    # ------------------------------------------------------------------ caveats
    story += [P("What this does and does not show", h2)] + bullets([
        "<b>It is a lead, not proof.</b> Every identification carries a small chance of being wrong, and a position supported by "
        "a single peptide is thin evidence. The peptides that matter should be checked by looking at their spectra.",
        "<b>One run, no repeats.</b> With a single run there is no way to see how repeatable the result is.",
        "<b>\"No peptide\" is not a negative result.</b> It only means the run did not happen to identify a peptide "
        "through that position.",
        "<b>The program counts peptides; it does not decide biology.</b> Other explanations (for example a different "
        "amino acid with the same mass, or an error in the genome annotation) are for the scientist to weigh."])
    story += [P("Sensible next steps", h2)] + bullets([
        "Open the spectra of the key peptides in the search software's viewer and check that the fragment pattern really "
        "contains the alternative letter.",
        "Repeat the experiment with replicates, so that the same peptides can be found more than once.",
        "To reach the uncovered positions, digest with a second enzyme (so that different peptides are produced) or run more material.",
        ("Compare the result with the standard-database search shown above, and look at the spectra at any position where the "
         "two searches disagree." if summ.get("has_b") else
         "Search the same data against the standard database alone and add that report, to see how the data behave when only "
         "the original reading is allowed.")])

    # ------------------------------------------------------------------ appendix
    story += [PageBreak(), P("Appendix A. All changed positions", h1)]
    rows = [[P("Protein", cellb), P("Pos.", cellb), P("Change", cellb), P("Verdict", cellb), P("Sequence around the position", cellb)]]
    for s in sites:
        rows.append([P(_esc(_short(s["ref_protein"])), cell), P(str(s["position"]), cell),
                     P(f"{_esc(s['ref_aa'])} to {_esc(s['alt_aa'])}", cell), P(_verdict_word(s["status"]), cell),
                     P(f"{_esc(s['context_ref'])}<br/>{_esc(s['context_alt'])}", ParagraphStyle(
                         "mono", parent=cell, fontName="Courier", fontSize=8))])
    story += [table(rows, [1.35 * inch, 0.5 * inch, 0.7 * inch, 1.2 * inch, 3.15 * inch]),
              P("The bracketed letter is the changed position: the first line is the standard sequence, the second the alternative.", small)]

    story += [P("Appendix B. Glossary", h1)]
    gloss = [("Peptide", "A short piece of a protein, made by cutting it with an enzyme. Mass spectrometry identifies peptides, "
                         "then infers the proteins."),
             ("Database (FASTA)", "The list of protein sequences the search software matches against."),
             ("Coverage", "How much was seen: either the share of proteins detected, or the share of a protein's letters that "
                          "fall inside identified peptides."),
             ("Position / site", "One letter (amino acid) in a protein, counted from the start."),
             ("Trypsin", "The usual enzyme; it cuts after K (lysine) and R (arginine), except before P (proline)."),
             ("Missed cleavage", "A cut site inside a peptide that the enzyme did not cut."),
             ("Intensity", "The measured signal strength. Useful for comparing peptides, not an absolute amount."),
             ("DIA-NN", "The software that identified and quantified the peptides from the raw data."),
             ("Fragment ion (b and y)", "When a peptide is broken in the instrument it splits into pieces. b ions hold the start of the "
                                        "sequence and y ions the end; their masses spell out the sequence."),
             ("MS/MS spectrum", "The list of fragment masses and strengths measured for one isolated peptide."),
             ("Retention time", "When a peptide comes off the chromatography column, in minutes after the run started."),
             ("Chance probability", "How likely it is that the observed number of fragment matches happened by accident. Smaller is more convincing.")]
    story += [table([[P(f"<b>{a}</b>", cell), P(b, cell)] for a, b in gloss], [1.5 * inch, 5.4 * inch], header=False)]

    story += [P("Appendix C. Methods and files", h1)]
    fl = [("Peptide report (main search)", files.get("report", "")), ("Reference FASTA", files.get("fasta_ref", "")),
          ("Alternative FASTA", files.get("fasta_alt", "") or "(variants were in the reference file)"),
          ("Standard-database search report", files.get("report_b", "") or "not provided")]
    story += [table([[P(a, cellb), P(_esc(b), cell)] for a, b in fl], [2.2 * inch, 4.7 * inch], header=False)]
    story += [Spacer(1, 6), P(
        f"Each protein in the alternative database was paired with its reference protein (by shared accession, or by matching ends) "
        f"and every position where the two differ was recorded ({summ['n_sites']} positions in {summ['n_proteins']} proteins; "
        f"{summ['n_unpaired']} alternative proteins could not be paired). Identified peptides (5 to 70 amino acids) were matched to "
        "positions: a peptide spans a position if the position lies within it. A peptide that also occurs in any other protein of the "
        "database is marked ambiguous and not counted. A position is judged \"alternative\" or \"standard\" if at least one unambiguous "
        "peptide of that version spans it and none of the other does, \"both\" if both versions have one, and \"no peptide\" otherwise. "
        "Coverage is computed against the reference proteins only, from the same peptide list. The identifications are as reliable as "
        "the search software's own false-discovery control; no further filtering was applied here.", body)]

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor(GREY))
        canvas.drawString(0.8 * inch, 0.5 * inch, f"Recode Detector {__version__} - database comparison report")
        canvas.drawRightString(letter[0] - 0.8 * inch, 0.5 * inch, f"Page {doc.page}")
        canvas.restoreState()

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter, leftMargin=0.8 * inch, rightMargin=0.8 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.8 * inch, title="Database comparison report",
                            author="Recode Detector")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()


def _verdict_word(status: str) -> str:
    return {"alternative": "Alternative version seen", "reference": "Standard version seen",
            "both": "Both seen", "none": "No peptide covers it"}[status]


def _short(pid: str) -> str:
    import re
    t = [x for x in pid.split("|") if x and not re.fullmatch(r"(CANONICAL|CGG2W|ALT|VARIANT)", x, flags=re.I)]
    return t[-1] if t else pid


def _median_len(cov: dict) -> str:
    items = sorted((int(k), v) for k, v in cov["pep_lengths"].items())
    total, run = sum(v for _, v in items), 0
    for k, v in items:
        run += v
        if run >= total / 2:
            return str(k)
    return "-"


def _spectra_section(spectra, sites, P, h2, cell, cellb, small, table, captioned, colors, inch, Drawing, Line, String, Rect, KeepTogether):
    res = spectra["results"]
    src = spectra.get("source", {})
    out = [P("Looking at the raw spectra", h2),
           P(f"The identifications above come from the search software. As an independent check, the program went back to the raw data "
             f"({_esc(src.get('name', 'the raw file'))}{', ' + _esc(src.get('instrument', '')) if src.get('instrument') else ''}) "
             f"and looked at the actual fragmentation spectra for the key peptides. For each peptide it finds the intact ion in the "
             f"survey scans, takes the fragmentation scans taken at the moment it elutes, and asks how many of the expected fragment "
             f"ions are really there. It then does the same for the <i>counterpart</i> peptide carrying the other letter at the changed "
             f"position. Fragments that contain the changed letter are the ones that tell the versions apart, so they carry the verdict. "
             f"The chance column is the probability that this many fragments would match by accident, corrected for picking the best of "
             f"many scans; smaller means more convincing. As a negative control, 20 shuffled decoy peptides of the same mass are scored on the "
             f"same scans; the Decoys column shows how many matched as well (it should be 0).")]
    word = {"supported": "Supported", "weak": "Weak", "not seen": "Not seen", "no MS1 signal": "No signal", "no MS/MS": "No MS/MS", "error": "Error"}
    rows = [[P("Protein", cellb), P("Pos", cellb), P("Version", cellb), P("Peptide", cellb), P("Elutes", cellb),
             P("Ions (changed)", cellb), P("Chance", cellb), P("Decoys", cellb), P("Verdict", cellb)]]
    for r in res:
        site = sites[r["site"]] if r.get("site") is not None and r["site"] < len(sites) else None
        m1, m2 = r.get("ms1") or {}, r.get("ms2") or {}
        bm = m2.get("best") or {}
        rows.append([P(_esc(_short(site["ref_protein"])) if site else "", cell), P(str(site["position"]) if site else "", cell),
                     P(_esc(r.get("version", "")), cell), P(_esc(r.get("modified", "")), cell),
                     P(f"{m1['apex_rt']:.2f}" if m1.get("apex_rt") else "-", cell),
                     P(f"{bm.get('n_matched', '-')} ({bm.get('n_disc', '-')})", cell),
                     P(("<1e-12" if bm["p_chance"] <= 1e-12 else f"{bm['p_chance']:.0e}") if bm.get("p_chance") is not None else "-", cell),
                     P(f"{bm['decoys']['as_good']}/{bm['decoys']['n']}" if bm.get("decoys") else "-", cell),
                     P(word.get(r.get("verdict"), r.get("verdict", "")), cell)])
    out += [table(rows, [0.95 * inch, 0.4 * inch, 0.85 * inch, 1.35 * inch, 0.55 * inch, 0.7 * inch, 0.6 * inch, 0.65 * inch, 0.75 * inch])]
    out += [P("Verdicts are heuristic: Supported needs at least two fragment ions that contain the changed letter, with a chance of "
              "1 in 1,000 or better, plus a good overall match and at most one decoy matching as well. Co-elution of the fragments (they should rise and fall together) is "
              "reported in the spreadsheet export. Because this is a data-dependent read of one run, treat it as a lead and confirm in "
              "the search software or a dedicated viewer.", small)]
    # mini annotated spectra for the strongest supported peptides
    good = [r for r in res if r.get("verdict") == "supported" and (r.get("ms2") or {}).get("best")]
    good.sort(key=lambda r: r["ms2"]["best"]["p_chance"])
    for r in good[:3]:
        bm = r["ms2"]["best"]
        w, h = 470, 150
        d = Drawing(w, h + 24)
        mz, it = bm["mz"], bm["int"]
        lo, hi = min(mz) - 5, max(mz) + 5
        mx = max(it) or 1
        d.add(Rect(30, 18, w - 40, h - 8, strokeColor=colors.HexColor("#cbd5e1"), fillColor=None, strokeWidth=.5))
        matched = {x["index"]: x for x in bm["matches"]}
        for i, (a, b) in enumerate(zip(mz, it)):
            x = 30 + (a - lo) / (hi - lo) * (w - 40)
            y = 18 + b / mx * (h - 30)
            hit = i in matched
            d.add(Line(x, 18, x, y, strokeColor=colors.HexColor("#b8b2a3" if not hit else (AMBER if matched[i]["disc"] else BLUE)), strokeWidth=1.3 if hit else .6))
            if hit:
                d.add(String(x - 6, y + 2, matched[i]["label"], fontSize=6.5, fontName="Helvetica", fillColor=colors.HexColor(GREY)))
        d.add(String(30, h + 12, f"{_esc(r['modified'])}  ({r['version']}), charge {r['charge']}+, scan {bm['scan']}, {r['ms1']['apex_rt']:.2f} min", fontSize=8, fontName="Helvetica-Bold"))
        out += [KeepTogether(captioned(d, "Matched fragment ions are labelled. Red = contains the changed residue; green = other matched ions; grey = unmatched."))]
    return out
